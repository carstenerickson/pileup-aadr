"""Tests for shard.py — Stage 3 fan-out shard manifest + merge."""
from __future__ import annotations

from pathlib import Path

import pytest

from pileup_aadr.errors import PileupAadrInternalError
from pileup_aadr.shard import (
    ShardSpec,
    build_shard_manifest,
    derive_shard_seed,
    merge_shard_eigenstrat,
)


def _write_sites(
    snp_path: Path,
    bed_path: Path,
    rows: list[tuple[str, str, int]],
) -> None:
    """Write minimal .snp + BED for testing.

    rows: (rsid, chrom_chr, pos) — chrom_chr is chr-prefixed (e.g. "chr1").
    """
    snp_path.parent.mkdir(parents=True, exist_ok=True)
    with open(snp_path, "w") as snp_out, open(bed_path, "w") as bed_out:
        for rsid, chrom_chr, pos in rows:
            chrom_int = chrom_chr[3:] if chrom_chr.startswith("chr") else chrom_chr
            snp_out.write(f"{rsid}\t{chrom_int}\t0.0\t{pos}\tA\tG\n")
            bed_out.write(f"{chrom_chr}\t{pos - 1}\t{pos}\n")


def _make_shard_output(prefix: Path, rsids: list[str], sample: str = "S") -> ShardSpec:
    """Write a fake pileupCaller shard triplet and return a ShardSpec pointing to it."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(f"{prefix}.geno", "w") as gh, open(f"{prefix}.snp", "w") as sh:
        for i, rsid in enumerate(rsids):
            gh.write("0\n")
            sh.write(f"{rsid}\t1\t0.0\t1000\tA\tG\n")
    Path(f"{prefix}.ind").write_text(f"{sample}\tU\tPop\n")
    return prefix


# --- build_shard_manifest ---


def test_manifest_groups_rows_by_chromosome(tmp_path: Path) -> None:
    """build_shard_manifest partitions sites into per-chromosome ShardSpecs."""
    snp = tmp_path / "sites.snp"
    bed = tmp_path / "sites.bed"
    _write_sites(snp, bed, [
        ("rs1", "chr1", 1000),
        ("rs2", "chr1", 2000),
        ("rs3", "chr22", 3000),
    ])
    manifest = build_shard_manifest(snp, bed, tmp_path / "shards", master_seed=42)

    assert len(manifest) == 2
    assert manifest[0].chromosome == "chr1"
    assert manifest[1].chromosome == "chr22"
    assert len(manifest[0].snp_path.read_text().splitlines()) == 2
    assert len(manifest[1].snp_path.read_text().splitlines()) == 1


def test_manifest_skips_chromosomes_with_no_sites(tmp_path: Path) -> None:
    """Chromosomes absent from the input are omitted from the manifest."""
    snp = tmp_path / "sites.snp"
    bed = tmp_path / "sites.bed"
    _write_sites(snp, bed, [("rs_x", "chrX", 5000)])
    manifest = build_shard_manifest(snp, bed, tmp_path / "shards", master_seed=1)

    assert len(manifest) == 1
    assert manifest[0].chromosome == "chrX"


def test_manifest_rejects_unknown_chromosome(tmp_path: Path) -> None:
    """BED row with chromosome not in CHROM_ORDER → PileupAadrInternalError."""
    snp_path = tmp_path / "sites.snp"
    bed_path = tmp_path / "sites.bed"
    snp_path.parent.mkdir(parents=True, exist_ok=True)
    snp_path.write_text("rs1\t1\t0.0\t1000\tA\tG\n")
    bed_path.write_text("chr1_random_alt\t999\t1000\n")

    with pytest.raises(PileupAadrInternalError, match="not in CHROM_ORDER"):
        build_shard_manifest(snp_path, bed_path, tmp_path / "shards", master_seed=42)


# --- derive_shard_seed ---


def test_seed_no_collision_across_adjacent_master_seeds() -> None:
    """Adjacent master seeds produce disjoint per-shard seed ranges across all shards."""
    seeds_m42 = {derive_shard_seed(42, i) for i in range(25)}
    seeds_m43 = {derive_shard_seed(43, i) for i in range(25)}
    assert len(seeds_m42 & seeds_m43) == 0
    assert derive_shard_seed(42, 1) == 42 * 1009 + 1


# --- merge_shard_eigenstrat ---


def test_merge_preserves_manifest_chrom_order(tmp_path: Path) -> None:
    """merge_shard_eigenstrat concatenates shards in manifest (CHROM_ORDER) order."""
    chr1_prefix = tmp_path / "chr1" / "call"
    chr22_prefix = tmp_path / "chr22" / "call"
    _make_shard_output(chr1_prefix, ["rs1", "rs2"])
    _make_shard_output(chr22_prefix, ["rs22"])

    manifest = [
        ShardSpec(
            shard_index=0, chromosome="chr1",
            bed_path=tmp_path / "chr1" / "s.bed",
            snp_path=tmp_path / "chr1" / "s.snp",
            output_prefix=chr1_prefix, seed=0,
        ),
        ShardSpec(
            shard_index=21, chromosome="chr22",
            bed_path=tmp_path / "chr22" / "s.bed",
            snp_path=tmp_path / "chr22" / "s.snp",
            output_prefix=chr22_prefix, seed=21,
        ),
    ]

    geno_path, snp_path, ind_path = merge_shard_eigenstrat(manifest, tmp_path / "merged")

    geno_lines = geno_path.read_text().splitlines()
    snp_lines = snp_path.read_text().splitlines()
    assert len(geno_lines) == 3
    assert snp_lines[0].startswith("rs1")
    assert snp_lines[2].startswith("rs22")


def test_merge_strict_raises_on_mismatched_line_counts(tmp_path: Path) -> None:
    """Shard with mismatched .geno/.snp row counts → ValueError (zip strict)."""
    prefix = tmp_path / "chr1" / "call"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{prefix}.geno").write_text("0\n0\n")  # 2 lines
    Path(f"{prefix}.snp").write_text("rs1\t1\t0.0\t1000\tA\tG\n")  # 1 line
    Path(f"{prefix}.ind").write_text("S\tU\tP\n")

    spec = ShardSpec(
        shard_index=0, chromosome="chr1",
        bed_path=prefix.parent / "b.bed",
        snp_path=prefix.parent / "s.snp",
        output_prefix=prefix, seed=0,
    )
    with pytest.raises(ValueError):
        merge_shard_eigenstrat([spec], tmp_path / "merged")
