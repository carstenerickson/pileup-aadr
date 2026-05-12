"""HG19_CHROM_LENGTHS / HG38_CHROM_LENGTHS integrity tests.

Maps to LLD test #30: catches a future maintainer who accidentally bumps a number
in either chromosome-length table. Values are byte-stable facts of the assembly
versions and should never change.
"""
from __future__ import annotations

from pileup_aadr.chrom_lengths import (
    CHROM_ORDER,
    HG19_CHROM_LENGTHS,
    HG38_CHROM_LENGTHS,
)


def test_hg19_table_has_25_entries() -> None:
    """22 autosomes + chrX + chrY + chrM."""
    assert len(HG19_CHROM_LENGTHS) == 25


def test_hg38_table_has_25_entries() -> None:
    assert len(HG38_CHROM_LENGTHS) == 25


def test_hg19_chr1_canonical() -> None:
    """The most-checked sentinel value across the codebase (build detection ±1Mb)."""
    assert HG19_CHROM_LENGTHS["chr1"] == 249_250_621


def test_hg38_chr1_canonical() -> None:
    assert HG38_CHROM_LENGTHS["chr1"] == 248_956_422


def test_hg19_chrM_canonical() -> None:
    """hg19 mtDNA is 16,571 bp (NC_001807, 'rCRS-clone')."""
    assert HG19_CHROM_LENGTHS["chrM"] == 16_571


def test_hg38_chrM_canonical() -> None:
    """hg38 mtDNA is 16,569 bp (NC_012920, 'rCRS') — distinct from hg19 by 2 bp."""
    assert HG38_CHROM_LENGTHS["chrM"] == 16_569


def test_chrom_order_is_22_autosomes_then_xym() -> None:
    """CHROM_ORDER drives ##contig emission order; must be numeric-then-XYM, not lex."""
    assert CHROM_ORDER[:22] == [f"chr{i}" for i in range(1, 23)]
    assert CHROM_ORDER[22:] == ["chrX", "chrY", "chrM"]


def test_chrom_order_covers_all_table_keys() -> None:
    """Every chrom in the length tables appears in CHROM_ORDER (no orphans)."""
    assert set(HG19_CHROM_LENGTHS) == set(CHROM_ORDER)
    assert set(HG38_CHROM_LENGTHS) == set(CHROM_ORDER)


def test_hg19_and_hg38_disagree_on_most_chroms() -> None:
    """Sanity: the two assemblies should differ in length on most chroms (build-detect basis)."""
    differing = sum(
        1
        for c in CHROM_ORDER
        if HG19_CHROM_LENGTHS[c] != HG38_CHROM_LENGTHS[c]
    )
    # All 22 autosomes + chrX + chrY + chrM differ between hg19 and hg38
    assert differing == 25, f"expected all 25 chroms to differ, got {differing}"
