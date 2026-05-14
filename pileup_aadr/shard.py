"""Chromosome sharding helpers for Stage 3 fan-out (v0.3).

Stage 3's `samtools mpileup | pileupCaller` chain is single-threaded by both
tools' nature. v0.3 parallelises Stage 3 by sharding the Stage 2 sites BED + .snp
by chromosome and running N parallel chains, then concatenating the per-shard
EIGENSTRAT triplets in CHROM_ORDER before Stage 4.

This module owns the shard manifest data structure and helpers. The orchestrator
owns the ThreadPoolExecutor and shard-result aggregation. Imports from
`chrom_lengths` (not `format_detect`) for the canonical-chrom map; `chrom_lengths`
is the source of truth for CHROM_ORDER, which drives Stage 2's BED/SNP ordering.
"""
from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from pathlib import Path

from .chrom_lengths import CHROM_ORDER
from .errors import PileupAadrInternalError

log = logging.getLogger(__name__)

# Prime multiplier for per-shard seed derivation. Larger than len(CHROM_ORDER)=25,
# so adjacent master seeds produce disjoint shard-seed ranges and cannot collide
# across shards (master=42 shard=1 → 42*1009+1=42379; master=43 shard=0 →
# 43*1009=43387; no overlap). Documented in JSON report config block.
_SHARD_SEED_MULTIPLIER: int = 1009


@dataclass(frozen=True)
class ShardSpec:
    """One Stage 3 shard — 1:1 with a chromosome in the per-chrom layout.

    The threads=1 path in `run_pileup_call_shards` does NOT use this class;
    it short-circuits directly to `run_pileup_call`.

    Attributes:
        shard_index: position in CHROM_ORDER for this shard's chromosome.
            Used as the `i` term in `seed = master * 1009 + i`.
        chromosome: canonical chr-prefixed name ("chr1".."chr22", "chrX", "chrY").
            "chrM" included only if the AADR panel has chrM sites (1240k does not).
        bed_path: shard-scoped BED (Stage 2 BED filtered to this chromosome only).
        snp_path: shard-scoped pileupCaller .snp (same filter; row order preserved).
        output_prefix: pileupCaller writes <output_prefix>.{geno,snp,ind} per shard.
        seed: pileupCaller --seed for this shard.
    """

    shard_index: int
    chromosome: str
    bed_path: Path
    snp_path: Path
    output_prefix: Path
    seed: int


def derive_shard_seed(master_seed: int, shard_index: int) -> int:
    """Per-shard seed: master_seed * 1009 + shard_index.

    The prime multiplier prevents shard-seed collisions across adjacent master
    seeds. For any fixed shard_index, seeds for consecutive master values are
    1009 apart; for any fixed master, seeds for consecutive shards are 1 apart.
    """
    return master_seed * _SHARD_SEED_MULTIPLIER + shard_index


def build_shard_manifest(
    sites_snp_path: Path,
    sites_bed_path: Path,
    shard_dir: Path,
    master_seed: int,
) -> list[ShardSpec]:
    """Partition the Stage 2 sites .snp + BED into per-chromosome shards.

    Stage 2 emits both files in CHROM_ORDER x pos_bp order (build_sites_vcf
    sorts pre-Picard; Picard preserves order; transform.py iterates the lifted
    VCF in file order), so grouping rows by chromosome in a single linear pass
    preserves per-chromosome position order. Empty chromosomes (no rows in input)
    are silently omitted from the manifest.

    Args:
        sites_snp_path: Stage 2 output .snp (numeric chrom in col 2; AADR-ordered)
        sites_bed_path: Stage 2 output BED (chr-prefixed in col 0; same row count)
        shard_dir: parent directory for per-shard files (e.g. <tempdir>/shards/)
        master_seed: pileupCaller master seed (from --seed CLI flag)

    Returns:
        list[ShardSpec], sorted by shard_index ascending. Non-empty (the
        orchestrator should have failed earlier if Stage 2 produced zero rows,
        but we guard with PileupAadrInternalError just in case).

    Raises:
        PileupAadrInternalError: a BED row has a chromosome not in CHROM_ORDER
            (Stage 2's alt_contig_filter should have removed it; defensive).
    """
    chrom_order_set = set(CHROM_ORDER)
    snp_rows_by_chrom: dict[str, list[str]] = {}
    bed_rows_by_chrom: dict[str, list[str]] = {}

    with open(sites_snp_path) as snp_in, open(sites_bed_path) as bed_in:
        for line_no, (snp_line, bed_line) in enumerate(
            zip(snp_in, bed_in, strict=True), start=1
        ):
            bed_chrom = bed_line.split("\t", 1)[0]
            if bed_chrom not in chrom_order_set:
                raise PileupAadrInternalError(
                    what=f"build_shard_manifest line {line_no}",
                    why=(
                        f"chromosome {bed_chrom!r} not in CHROM_ORDER; "
                        "Stage 2 alt_contig_filter should have removed this"
                    ),
                    fix=(
                        "Verify Stage 2's _CANONICAL_CHROM_RE; report a bug if "
                        "the chrom is canonical but absent from CHROM_ORDER"
                    ),
                )
            snp_rows_by_chrom.setdefault(bed_chrom, []).append(snp_line)
            bed_rows_by_chrom.setdefault(bed_chrom, []).append(bed_line)

    manifest: list[ShardSpec] = []
    for shard_index, chromosome in enumerate(CHROM_ORDER):
        if chromosome not in bed_rows_by_chrom:
            log.debug("shard %d (%s): no sites; skipped", shard_index, chromosome)
            continue
        chrom_dir = shard_dir / chromosome
        chrom_dir.mkdir(parents=True, exist_ok=True)
        snp_path = chrom_dir / "sites.snp"
        bed_path = chrom_dir / "sites.bed"
        with open(snp_path, "w") as snp_out:
            snp_out.writelines(snp_rows_by_chrom[chromosome])
        with open(bed_path, "w") as bed_out:
            bed_out.writelines(bed_rows_by_chrom[chromosome])
        manifest.append(ShardSpec(
            shard_index=shard_index,
            chromosome=chromosome,
            bed_path=bed_path,
            snp_path=snp_path,
            output_prefix=chrom_dir / "call",
            seed=derive_shard_seed(master_seed, shard_index),
        ))

    if not manifest:
        raise PileupAadrInternalError(
            what="build_shard_manifest",
            why=(
                "no shards produced (Stage 2 output was empty or all rows had "
                "non-CHROM_ORDER chromosomes)"
            ),
            fix="Inspect Stage 2 counters; this should not fire if Stage 3 was reached",
        )

    log.info(
        "Built shard manifest: %d shards (%s)",
        len(manifest),
        ", ".join(s.chromosome for s in manifest),
    )
    return manifest


def merge_shard_eigenstrat(
    manifest: list[ShardSpec],
    merged_prefix: Path,
) -> tuple[Path, Path, Path]:
    """Concatenate per-shard EIGENSTRAT triplets into a single merged set.

    Concatenates .geno and .snp in manifest order (= CHROM_ORDER). All shards
    share the same --sample-name / --sex / --pop-name CLI args, so their .ind
    files are byte-identical; we copy shard 0's .ind after a defensive equality
    check.

    The threads=1 short-circuit in run_pileup_call_shards never calls this
    function — that path writes directly to merged_prefix.

    Args:
        manifest: shards in CHROM_ORDER (output of build_shard_manifest).
        merged_prefix: where to write the merged .geno + .snp + .ind.

    Returns:
        (merged_geno_path, merged_snp_path, merged_ind_path)

    Raises:
        ValueError: a shard's .geno and .snp have mismatched line counts
            (from zip strict=True; pileupCaller should always emit aligned files).
        PileupAadrInternalError: shard .ind files are not byte-identical
            (indicates inconsistent CLI args passed to shards).
    """
    merged_geno_path = Path(f"{merged_prefix}.geno")
    merged_snp_path = Path(f"{merged_prefix}.snp")
    merged_ind_path = Path(f"{merged_prefix}.ind")
    merged_geno_path.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    with open(merged_geno_path, "w") as out_geno, open(merged_snp_path, "w") as out_snp:
        for shard in manifest:
            shard_geno = Path(f"{shard.output_prefix}.geno")
            shard_snp = Path(f"{shard.output_prefix}.snp")
            shard_rows = 0
            with open(shard_geno) as in_geno, open(shard_snp) as in_snp:
                for geno_line, snp_line in zip(in_geno, in_snp, strict=True):
                    out_geno.write(geno_line)
                    out_snp.write(snp_line)
                    shard_rows += 1
            log.debug("merged shard %s: %d rows", shard.chromosome, shard_rows)
            total_rows += shard_rows

    shard0_ind_bytes = Path(f"{manifest[0].output_prefix}.ind").read_bytes()
    for shard in manifest[1:]:
        shard_ind_bytes = Path(f"{shard.output_prefix}.ind").read_bytes()
        if shard_ind_bytes != shard0_ind_bytes:
            raise PileupAadrInternalError(
                what=f"merge_shard_eigenstrat .ind mismatch ({shard.chromosome})",
                why=(
                    f"shard {shard.chromosome}'s .ind differs from shard "
                    f"{manifest[0].chromosome}'s — all shards must have the "
                    "same --sample-name / --pop-name / --sex"
                ),
                fix=(
                    "Inspect per-shard .ind files; indicates the orchestrator "
                    "passed inconsistent args to shards"
                ),
            )
    shutil.copy(Path(f"{manifest[0].output_prefix}.ind"), merged_ind_path)

    log.info(
        "Merged %d shards into %s: %d total rows",
        len(manifest), merged_prefix, total_rows,
    )
    return merged_geno_path, merged_snp_path, merged_ind_path


__all__ = [
    "ShardSpec",
    "build_shard_manifest",
    "derive_shard_seed",
    "merge_shard_eigenstrat",
]
