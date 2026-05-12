"""Stage 2: convert Picard's lifted VCF to pileupCaller .snp + mpileup BED.

Two artifacts written to `<tempdir>/transform/`:

1. **pileupCaller .snp** — 6-col EIGENSOFT (numeric-chrom, AADR-style)
2. **mpileup BED** — 3-col 0-based BED (chr-prefixed, matches modern hg38 BAM @SQ)

The alt-contig filter is load-bearing: pileupCaller's Haskell parser raises an
uncatchable SeqFormatException on alt/decoy contigs 5-10 minutes into Stage 3,
so we drop them here per HLD §"Alt-contig filter".
"""
from __future__ import annotations

import logging
import re
import time
from pathlib import Path
from typing import Final

import pysam

from .counters import Stage2TransformCounters

log = logging.getLogger(__name__)

# Canonical-chrom regex: accepts only chroms pileupCaller's parseSnpFile handles
# (verified against pileupCaller.hs source). Everything else is dropped.
_CANONICAL_CHROM_RE: Final[re.Pattern[str]] = re.compile(
    r"^(chr)?([0-9]{1,2}|X|Y|MT|M)$"
)

# AADR/EIGENSOFT chromosome encoding for the .snp file's column 2 (numeric form).
# pileupCaller auto-converts internally; we emit numeric to match AADR convention.
_CHROM_TO_NUMERIC: Final[dict[str, str]] = {
    **{f"chr{i}": str(i) for i in range(1, 23)},
    "chrX": "23",
    "chrY": "24",
    "chrM": "90",
}


def build_pileupcaller_snp_and_bed(
    lifted_vcf_path: Path,
    output_snp_path: Path,
    output_bed_path: Path,
    alt_contig_filter: bool = True,
) -> Stage2TransformCounters:
    """Convert Picard's lifted VCF to pileupCaller .snp + mpileup BED.

    Args:
        lifted_vcf_path: from `lift.lift_aadr_sites` (Picard's OUTPUT)
        output_snp_path: where to write the .snp file
            (typically `<tempdir>/transform/aadr_hg38.snp`)
        output_bed_path: where to write the BED file
            (typically `<tempdir>/transform/aadr_hg38.bed`)
        alt_contig_filter: drop alt-haplotype + decoy contigs (default True; almost
            always wanted — pileupCaller crashes mid-Stage-3 without this)

    Returns:
        Stage2TransformCounters with alt_contig_drops + output_sites + wallclock.
    """
    t0 = time.perf_counter()
    alt_contig_drops = 0
    output_sites = 0

    output_snp_path.parent.mkdir(parents=True, exist_ok=True)
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)

    with (
        pysam.VariantFile(str(lifted_vcf_path)) as vcf,
        open(output_snp_path, "w") as snp_out,
        open(output_bed_path, "w") as bed_out,
    ):
        for rec in vcf:
            chrom = rec.chrom
            if alt_contig_filter and _CANONICAL_CHROM_RE.match(chrom) is None:
                alt_contig_drops += 1
                log.debug("Dropping alt-contig site: %s:%d", chrom, rec.pos)
                continue
            # Defensive: Picard already emits chr-prefixed but normalize anyway
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            chrom_numeric = _CHROM_TO_NUMERIC.get(chrom)
            if chrom_numeric is None:
                # Reachable only with alt_contig_filter=False on a non-canonical chrom
                # that happened to slip the regex (defensive double-check)
                alt_contig_drops += 1
                continue
            aadr_rs = rec.info.get("AADR_RS")
            rsid = aadr_rs if aadr_rs else rec.id  # fallback to ID col if INFO missing
            ref = rec.ref
            alts = rec.alts or ()
            if len(alts) != 1:
                # Picard never emits multi-allelic from a biallelic input; defensive
                log.warning(
                    "Skipping multi-allelic lifted record at %s:%d (%d ALTs)",
                    chrom, rec.pos, len(alts),
                )
                continue
            alt = alts[0]
            # .snp row: rsid \t chrom_numeric \t 0.0 \t pos_bp \t REF \t ALT
            # Genetic distance: 0.0 — Picard didn't carry Morgans through; Stage 4
            # rejoin uses AADR's Morgans verbatim, so this column is unused downstream.
            snp_out.write(
                f"{rsid}\t{chrom_numeric}\t0.0\t{rec.pos}\t{ref}\t{alt}\n"
            )
            # BED row: chrom (chrN) \t pos-1 (0-based) \t pos (1-based end)
            bed_out.write(f"{chrom}\t{rec.pos - 1}\t{rec.pos}\n")
            output_sites += 1

    wallclock = time.perf_counter() - t0
    log.info(
        "Stage 2 complete: %d sites written to .snp + .bed; "
        "%d alt-contig sites dropped; wallclock %.1fs",
        output_sites, alt_contig_drops, wallclock,
    )
    return Stage2TransformCounters(
        wallclock_seconds=wallclock,
        alt_contig_drops=alt_contig_drops,
        output_sites=output_sites,
    )


__all__ = ["build_pileupcaller_snp_and_bed"]
