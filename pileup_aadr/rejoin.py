"""Stage 4: rsID rejoin + SwappedAlleles dosage inversion.

The most algorithmically dense stage. Walks pileupCaller's per-variant output,
looks up the corresponding AADR row + lifted INFO record, applies SWAP_DOSAGE
inversion when Picard flagged a swap, and writes the final EIGENSTRAT triplet
in AADR's hg19 frame. Also computes the coverage gate's per-chrom breakdown.

For the no-lift fast path (AADR build == BAM build), Stage 4 is trivialized:
pileupCaller's output IS the final triplet (no rejoin needed). The
`_no_lift_fast_path_finalize` helper just copies the triplet, walks .geno+.snp
once for counters, overrides .ind with the user --sex, and constructs the
sidecar.
"""
from __future__ import annotations

import logging
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import pandas as pd
import pysam

from .counters import CoverageCounters, Stage4RejoinCounters
from .errors import PileupAadrInternalError
from .format_detect import normalize_chrom

log = logging.getLogger(__name__)

# .geno on-disk encoding (single ASCII chars, one per line per sample per variant)
GENO_HOM_REF: Final[str] = "0"
GENO_HET: Final[str] = "1"
GENO_HOM_ALT: Final[str] = "2"
GENO_MISSING: Final[str] = "9"

# Dosage inversion table for Picard's SwappedAlleles. Keyed by str char to match
# .geno encoding directly (no int lookups, no parse-and-reformat).
SWAP_DOSAGE: Final[dict[str, str]] = {"0": "2", "1": "1", "2": "0", "9": "9"}

# Coverage gate denominator chroms — autosomes only. Sex chroms appear in the
# per-chrom breakdown but don't enter the gate fraction.
_AUTOSOME_CHROMS: Final[frozenset[str]] = frozenset(f"chr{i}" for i in range(1, 23))

_EIG_SNP_COLS: Final[int] = 6

# Flat lookup dict field layout for Stage 4 hot loop (v0.3 item D).
# {rsid: (chrom_int, gen_morgans, pos_bp, ref, alt, is_swapped)}
AADR_LOOKUP_FIELDS: Final[tuple[str, ...]] = (
    "chrom_int", "gen_morgans", "pos_bp", "ref", "alt", "is_swapped"
)
AadrLookupValue = tuple[str, float, int, str, str, bool]
AadrLookup = dict[str, AadrLookupValue]


@dataclass
class RejoinOutput:
    """Bundle of artifacts produced by Stage 4 (or the no-lift finalizer)."""

    stage_4_counters: Stage4RejoinCounters
    coverage_counters: CoverageCounters
    per_variant_rows: list[dict[str, Any]] = field(default_factory=list)
    pseudohaploid_sidecar: dict[str, Any] = field(default_factory=dict)


def build_swap_lookup(lifted_vcf_path: Path) -> dict[str, bool]:
    """Read Picard's lifted VCF and build {rsid: bool} for SwappedAlleles flag.

    rsid comes from the AADR_RS INFO field (preserved through Picard's lift).
    Called by the orchestrator after Stage 1; the result feeds build_merged_lookup
    before Stage 4.
    """
    swap_map: dict[str, bool] = {}
    with pysam.VariantFile(str(lifted_vcf_path)) as vcf:
        for rec in vcf:
            aadr_rs = rec.info.get("AADR_RS")
            if not aadr_rs:
                # Defensive: AADR_RS should always be present after Stage 1
                continue
            # SwappedAlleles is a Flag-type INFO field (no value, just present/absent)
            swap_map[aadr_rs] = "SwappedAlleles" in rec.info
    return swap_map


def build_merged_lookup(
    aadr_df: pd.DataFrame,
    swap_lookup: dict[str, bool],
) -> AadrLookup:
    """Build flat {rsid: (chrom_int, gen_morgans, pos_bp, ref, alt, is_swapped)} dict.

    Single itertuples pass over the 1.2M-row DataFrame; swap flags merged inline
    from swap_lookup. Eliminates pandas .loc overhead and the separate swap_lookup
    query from the Stage 4 hot loop.

    Args:
        aadr_df: AADR DataFrame with index=rsid.
        swap_lookup: {rsid: bool} from build_swap_lookup (empty dict for no-lift path).

    Returns:
        AadrLookup ready for Stage 4's rejoin_aadr_frame.
    """
    lookup: AadrLookup = {}
    for row in aadr_df.itertuples():
        rsid = row.Index
        lookup[rsid] = (
            row.chrom_int,
            row.gen_morgans,
            row.pos_bp,
            row.ref,
            row.alt,
            swap_lookup.get(rsid, False),
        )
    return lookup


def rejoin_aadr_frame(
    pileupcaller_eig_prefix: Path,
    aadr_lookup: AadrLookup,
    output_prefix: Path,
    sample_name: str,
    pop_name: str,
    *,
    sex: str = "U",
    emit_per_variant_rows: bool = False,
) -> RejoinOutput:
    """Rejoin pileupCaller output to AADR's hg19 frame; invert dosage at SwappedAlleles.

    Args:
        pileupcaller_eig_prefix: pileupCaller's <prefix> (reads .geno + .snp).
        aadr_lookup: flat lookup dict from build_merged_lookup (built by orchestrator
            after Stage 3 so swap_lookup Future is resolved before this call).
        output_prefix: where to write final `<prefix>.{geno,snp,ind}`.
        sample_name: IID for output .ind.
        pop_name: POP for output .ind.
        sex: SEX for output .ind ("M", "F", or "U"; default "U").
        emit_per_variant_rows: populate `RejoinOutput.per_variant_rows` for `--report-tsv`.

    Returns:
        RejoinOutput with counters + per-variant rows + sidecar dict.

    Raises:
        PileupAadrInternalError: malformed pileupCaller .snp row (should never
            fire — pileupCaller emits clean rows; defensive).
    """
    t0 = time.perf_counter()

    pc_geno_path = Path(f"{pileupcaller_eig_prefix}.geno")
    pc_snp_path = Path(f"{pileupcaller_eig_prefix}.snp")

    rsid_matched = 0
    ref_alt_swap_count = 0
    allele_mismatch_drops = 0
    output_variants = 0
    per_chrom_call_count: dict[str, int] = {}
    per_chrom_total: dict[str, int] = {}
    het_count = 0
    non_missing_autosomal = 0
    per_variant_rows: list[dict[str, Any]] = []

    out_geno_path = Path(f"{output_prefix}.geno")
    out_snp_path = Path(f"{output_prefix}.snp")
    out_geno_path.parent.mkdir(parents=True, exist_ok=True)

    with (
        open(pc_geno_path) as pc_geno,
        open(pc_snp_path) as pc_snp,
        open(out_geno_path, "w") as out_geno,
        open(out_snp_path, "w") as out_snp,
    ):
        for line_no, (geno_line, snp_line) in enumerate(
            zip(pc_geno, pc_snp, strict=True), start=1
        ):
            geno_char = geno_line.rstrip()
            snp_parts = snp_line.split()
            if len(snp_parts) != _EIG_SNP_COLS:
                raise PileupAadrInternalError(
                    what=f"pileupCaller .snp parse at line {line_no}",
                    why=f"expected 6 cols, got {len(snp_parts)}: {snp_line!r}",
                    fix="Verify pileupCaller output integrity; this should never fire",
                )
            rsid, _, _, _, lifted_ref, lifted_alt = snp_parts

            entry = aadr_lookup.get(rsid)
            if entry is None:
                # rsID in pileupCaller output but missing from AADR. Defensive:
                # we built the .snp from AADR so this shouldn't happen, but a
                # corrupt intermediate or upstream rename could trigger it.
                log.warning(
                    "rsid %s in pileupCaller output but missing from AADR; skipping", rsid
                )
                continue
            rsid_matched += 1
            aadr_chrom_int, aadr_gen, aadr_pos, aadr_ref, aadr_alt, was_swapped = entry

            if was_swapped:
                if not (lifted_ref == aadr_alt and lifted_alt == aadr_ref):
                    allele_mismatch_drops += 1
                    log.warning(
                        "rsid %s flagged SwappedAlleles but alleles don't match swap "
                        "pattern: lifted=%s/%s vs AADR=%s/%s",
                        rsid, lifted_ref, lifted_alt, aadr_ref, aadr_alt,
                    )
                    if emit_per_variant_rows:
                        per_variant_rows.append({
                            "aadr_id": rsid, "chrom_hg19": aadr_chrom_int,
                            "pos_hg19": aadr_pos, "ref_hg19": aadr_ref,
                            "alt_hg19": aadr_alt,
                            "action": "dropped_allele_mismatch",
                        })
                    continue
                new_geno = SWAP_DOSAGE[geno_char]
                ref_alt_swap_count += 1
                action = "swap"
            else:
                if not (lifted_ref == aadr_ref and lifted_alt == aadr_alt):
                    allele_mismatch_drops += 1
                    log.warning(
                        "rsid %s no swap flag but alleles differ: "
                        "lifted=%s/%s vs AADR=%s/%s",
                        rsid, lifted_ref, lifted_alt, aadr_ref, aadr_alt,
                    )
                    if emit_per_variant_rows:
                        per_variant_rows.append({
                            "aadr_id": rsid, "chrom_hg19": aadr_chrom_int,
                            "pos_hg19": aadr_pos, "ref_hg19": aadr_ref,
                            "alt_hg19": aadr_alt,
                            "action": "dropped_allele_mismatch",
                        })
                    continue
                new_geno = geno_char
                action = "passthrough" if new_geno != GENO_MISSING else "missing_call"

            out_geno.write(new_geno + "\n")
            out_snp.write(
                f"{rsid}\t{aadr_chrom_int}\t{aadr_gen}\t{aadr_pos}\t{aadr_ref}\t{aadr_alt}\n"
            )
            output_variants += 1

            chrom_chr = normalize_chrom(aadr_chrom_int) or aadr_chrom_int
            per_chrom_total[chrom_chr] = per_chrom_total.get(chrom_chr, 0) + 1
            if new_geno != GENO_MISSING:
                per_chrom_call_count[chrom_chr] = per_chrom_call_count.get(chrom_chr, 0) + 1
                if chrom_chr in _AUTOSOME_CHROMS:
                    non_missing_autosomal += 1
                if new_geno == GENO_HET:
                    het_count += 1

            if emit_per_variant_rows:
                per_variant_rows.append({
                    "aadr_id": rsid, "chrom_hg19": aadr_chrom_int,
                    "pos_hg19": aadr_pos, "ref_hg19": aadr_ref,
                    "alt_hg19": aadr_alt,
                    "action": action,
                })

    out_ind_path = Path(f"{output_prefix}.ind")
    with open(out_ind_path, "w") as out_ind:
        out_ind.write(f"{sample_name}\t{sex}\t{pop_name}\n")

    autosomal_aadr_count = sum(
        per_chrom_total.get(f"chr{i}", 0) for i in range(1, 23)
    )
    coverage_fraction = (
        non_missing_autosomal / autosomal_aadr_count if autosomal_aadr_count else 0.0
    )

    sidecar = _build_sidecar(
        sample_name=sample_name,
        het_count=het_count,
        non_missing_autosomal=non_missing_autosomal,
        no_lift=False,
    )

    wallclock = time.perf_counter() - t0
    log.info(
        "Stage 4 complete: %d output variants (%d swapped, %d dropped); wallclock %.1fs",
        output_variants, ref_alt_swap_count, allele_mismatch_drops, wallclock,
    )

    return RejoinOutput(
        stage_4_counters=Stage4RejoinCounters(
            wallclock_seconds=wallclock,
            rsid_matched=rsid_matched,
            ref_alt_swap_count=ref_alt_swap_count,
            allele_mismatch_drops=allele_mismatch_drops,
            output_variants=output_variants,
        ),
        coverage_counters=CoverageCounters(
            non_missing_autosomal_calls=non_missing_autosomal,
            coverage_fraction=round(coverage_fraction, 4),
            coverage_warning=False,  # set by orchestrator after gate evaluation
            per_chrom_call_count=per_chrom_call_count,
        ),
        per_variant_rows=per_variant_rows,
        pseudohaploid_sidecar=sidecar,
    )


def _no_lift_fast_path_finalize(
    pileupcaller_eig_prefix: Path,
    aadr_df: pd.DataFrame,
    output_prefix: Path,
    sample_name: str,
    pop_name: str,
    *,
    sex: str = "U",
    emit_per_variant_rows: bool = False,
) -> RejoinOutput:
    """Finalizer for the no-lift fast path (AADR build == BAM build).

    No SwappedAlleles handling, no rejoin — pileupCaller's output IS the AADR-
    frame triplet. Just copy + override .ind + compute coverage counters.

    Args:
        pileupcaller_eig_prefix: pileupCaller's <prefix> (reads .geno + .snp).
        aadr_df: AADR DataFrame (used for autosomal-row denominator).
        output_prefix: final EIGENSTRAT prefix.
        sample_name: IID for output .ind.
        pop_name: POP for output .ind.
        sex: SEX for output .ind ("M", "F", or "U"; default "U").
        emit_per_variant_rows: populate `RejoinOutput.per_variant_rows`.

    Returns:
        RejoinOutput with counters + per-variant rows + sidecar dict.
    """
    t0 = time.perf_counter()

    pc_geno_path = Path(f"{pileupcaller_eig_prefix}.geno")
    pc_snp_path = Path(f"{pileupcaller_eig_prefix}.snp")

    out_geno_path = Path(f"{output_prefix}.geno")
    out_snp_path = Path(f"{output_prefix}.snp")
    out_geno_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(pc_geno_path, out_geno_path)
    shutil.copy2(pc_snp_path, out_snp_path)

    output_variants = 0
    het_count = 0
    non_missing_autosomal = 0
    per_chrom_call_count: dict[str, int] = {}
    per_variant_rows: list[dict[str, Any]] = []

    with open(pc_geno_path) as gh, open(pc_snp_path) as sh:
        for geno_line, snp_line in zip(gh, sh, strict=True):
            geno_char = geno_line.rstrip()
            snp_parts = snp_line.split()
            if len(snp_parts) != _EIG_SNP_COLS:
                # Defensive — pileupCaller emits clean rows
                log.warning(
                    "malformed .snp row at variant %d: %r",
                    output_variants + 1, snp_line,
                )
                continue
            rsid, chrom_int, _, pos, ref, alt = snp_parts
            chrom_chr = normalize_chrom(chrom_int) or chrom_int
            output_variants += 1

            if geno_char != GENO_MISSING:
                per_chrom_call_count[chrom_chr] = per_chrom_call_count.get(chrom_chr, 0) + 1
                if chrom_chr in _AUTOSOME_CHROMS:
                    non_missing_autosomal += 1
                if geno_char == GENO_HET:
                    het_count += 1

            if emit_per_variant_rows:
                per_variant_rows.append({
                    "aadr_id": rsid, "chrom_hg19": chrom_int, "pos_hg19": int(pos),
                    "ref_hg19": ref, "alt_hg19": alt,
                    "action": "missing_call" if geno_char == GENO_MISSING else "passthrough",
                })

    # Override .ind with user-supplied sex (pileupCaller always writes SEX=U)
    out_ind_path = Path(f"{output_prefix}.ind")
    with open(out_ind_path, "w") as out_ind:
        out_ind.write(f"{sample_name}\t{sex}\t{pop_name}\n")

    sidecar = _build_sidecar(
        sample_name=sample_name,
        het_count=het_count,
        non_missing_autosomal=non_missing_autosomal,
        no_lift=True,
    )

    autosomal_aadr_count = int(
        aadr_df[aadr_df["chrom_int"].isin([str(i) for i in range(1, 23)])].shape[0]
    )
    coverage_fraction = (
        non_missing_autosomal / autosomal_aadr_count if autosomal_aadr_count else 0.0
    )

    wallclock = time.perf_counter() - t0
    log.info(
        "No-lift finalize complete: %d output variants (no swap handling); "
        "wallclock %.1fs",
        output_variants, wallclock,
    )

    return RejoinOutput(
        stage_4_counters=Stage4RejoinCounters(
            wallclock_seconds=wallclock,
            rsid_matched=output_variants,  # all rsIDs match by definition (no walk-and-lookup)
            ref_alt_swap_count=0,  # no swaps possible without lift
            allele_mismatch_drops=0,
            output_variants=output_variants,
        ),
        coverage_counters=CoverageCounters(
            non_missing_autosomal_calls=non_missing_autosomal,
            coverage_fraction=round(coverage_fraction, 4),
            coverage_warning=False,  # set by orchestrator after gate evaluation
            per_chrom_call_count=per_chrom_call_count,
        ),
        per_variant_rows=per_variant_rows,
        pseudohaploid_sidecar=sidecar,
    )


def _build_sidecar(
    *,
    sample_name: str,
    het_count: int,
    non_missing_autosomal: int,
    no_lift: bool,
) -> dict[str, Any]:
    """Build the PSEUDOHAPLOID sidecar dict (consumed by pgen-samplebind)."""
    note = (
        "no-lift fast path (AADR build == BAM build); single-BAM --randomDiploid "
        "output is pseudohaploid by construction"
        if no_lift
        else "single-BAM --randomDiploid output is pseudohaploid by construction"
    )
    return {
        "schema_version": 1,
        "samples": {
            sample_name: {
                "pseudohaploid": 1,
                "het_count": het_count,
                "non_missing_autosomal_count": non_missing_autosomal,
                "het_rate": (
                    het_count / non_missing_autosomal if non_missing_autosomal else 0.0
                ),
                "source": "pileup-aadr-extract",
                "calling_mode": "randomDiploid",
                "note": note,
            },
        },
    }


__all__ = [
    "AADR_LOOKUP_FIELDS",
    "AadrLookup",
    "AadrLookupValue",
    "GENO_HET",
    "GENO_HOM_ALT",
    "GENO_HOM_REF",
    "GENO_MISSING",
    "SWAP_DOSAGE",
    "RejoinOutput",
    "_no_lift_fast_path_finalize",
    "build_merged_lookup",
    "build_swap_lookup",
    "rejoin_aadr_frame",
]
