"""Tests for sites_vcf.py — pre-Stage 1 VCF construction from AADR DataFrame.

Maps to LLD test #29 (sites-VCF construction parses cleanly via bcftools/pysam,
##contig lines match table, AADR_RS INFO field present, chrom names use chrN form).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pileup_aadr.chrom_lengths import HG19_CHROM_LENGTHS
from pileup_aadr.format_detect import parse_aadr_snp
from pileup_aadr.sites_vcf import build_sites_vcf

# --- Basic construction ---


def test_build_sites_vcf_from_chr22_slice(aadr_chr22_slice: Path, tmp_path: Path) -> None:
    """The Day-1 AADR slice (50 chr22 sites) round-trips into a parseable VCF."""
    df = parse_aadr_snp(aadr_chr22_slice)
    out_vcf = tmp_path / "aadr_sites.vcf"
    counters = build_sites_vcf(df, out_vcf, aadr_build="hg19")

    assert out_vcf.exists()
    # The Day-1 slice has 1 palindrome (G>C) — verified in test_inspect_palindrome_count
    # via the fixture's allele distribution. So 50 input rows - 1 palindrome = 49 written.
    assert counters.rows_written == 49
    assert counters.palindrome_drops == 1
    assert counters.non_snp_drops == 0


def test_build_sites_vcf_emits_valid_vcf_v4_2_header(
    aadr_chr22_slice: Path, tmp_path: Path
) -> None:
    """Header has fileformat + ##contig + AADR_RS INFO + 8-col CHROM line."""
    df = parse_aadr_snp(aadr_chr22_slice)
    out_vcf = tmp_path / "out.vcf"
    build_sites_vcf(df, out_vcf, aadr_build="hg19")

    content = out_vcf.read_text()
    lines = content.splitlines()

    assert lines[0] == "##fileformat=VCFv4.2"
    assert any(line.startswith("##contig=<ID=chr22,length=") for line in lines)
    assert any("AADR_RS" in line and line.startswith("##INFO=") for line in lines)
    # CHROM header line
    chrom_line = next(line for line in lines if line.startswith("#CHROM"))
    assert chrom_line == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


def test_build_sites_vcf_chrom_chr22_length_matches_hg19_table(
    aadr_chr22_slice: Path, tmp_path: Path
) -> None:
    """##contig length for chr22 matches HG19_CHROM_LENGTHS exactly."""
    df = parse_aadr_snp(aadr_chr22_slice)
    out_vcf = tmp_path / "out.vcf"
    build_sites_vcf(df, out_vcf, aadr_build="hg19")

    content = out_vcf.read_text()
    expected = HG19_CHROM_LENGTHS["chr22"]
    assert f"##contig=<ID=chr22,length={expected}>" in content


def test_build_sites_vcf_emits_chrN_form_not_numeric(  # noqa: N802 — chrN refers to the contig naming convention
    aadr_chr22_slice: Path, tmp_path: Path
) -> None:
    """AADR's numeric chrom (22) must be normalized to chr-prefixed (chr22) on output."""
    df = parse_aadr_snp(aadr_chr22_slice)
    out_vcf = tmp_path / "out.vcf"
    build_sites_vcf(df, out_vcf, aadr_build="hg19")

    data_lines = [
        line for line in out_vcf.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    for line in data_lines:
        chrom = line.split("\t", 1)[0]
        assert chrom.startswith("chr"), f"row uses non-chrN form: {line!r}"


def test_build_sites_vcf_aadr_rs_info_per_row(
    aadr_chr22_slice: Path, tmp_path: Path
) -> None:
    """Every data row has AADR_RS=<rsid> in the INFO column for Stage-4 lookup."""
    df = parse_aadr_snp(aadr_chr22_slice)
    out_vcf = tmp_path / "out.vcf"
    build_sites_vcf(df, out_vcf, aadr_build="hg19")

    data_lines = [
        line for line in out_vcf.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    for line in data_lines:
        rsid = line.split("\t")[2]  # ID column
        info = line.split("\t")[7]  # INFO column
        assert info == f"AADR_RS={rsid}"


# --- Palindrome filter ---


def _df_from_rows(rows: list[tuple[str, str, float, int, str, str]]) -> pd.DataFrame:
    df = pd.DataFrame(
        rows,
        columns=["rsid", "chrom_int", "gen_morgans", "pos_bp", "ref", "alt"],
    )
    return df.set_index("rsid", verify_integrity=True)


def test_palindrome_filter_drops_at_default(tmp_path: Path) -> None:
    """A/T and C/G (and reverse) drop with default palindrome_filter=True."""
    df = _df_from_rows([
        ("rs_at", "22", 0.0, 1000, "A", "T"),
        ("rs_ta", "22", 0.0, 2000, "T", "A"),
        ("rs_cg", "22", 0.0, 3000, "C", "G"),
        ("rs_gc", "22", 0.0, 4000, "G", "C"),
        ("rs_ag", "22", 0.0, 5000, "A", "G"),  # not a palindrome
    ])
    counters = build_sites_vcf(df, tmp_path / "out.vcf", aadr_build="hg19")
    assert counters.palindrome_drops == 4
    assert counters.rows_written == 1


def test_palindrome_filter_opt_out(tmp_path: Path) -> None:
    """palindrome_filter=False passes A/T + C/G through (for diagnostic mode)."""
    df = _df_from_rows([
        ("rs_at", "22", 0.0, 1000, "A", "T"),
        ("rs_ag", "22", 0.0, 2000, "A", "G"),
    ])
    counters = build_sites_vcf(
        df, tmp_path / "out.vcf", aadr_build="hg19", palindrome_filter=False
    )
    assert counters.palindrome_drops == 0
    assert counters.rows_written == 2


# --- Non-SNP filter ---


def test_non_snp_filter_drops_indels(tmp_path: Path) -> None:
    """Multi-character REF/ALT (indels) drop with default non_snp_filter=True."""
    df = _df_from_rows([
        ("rs_indel1", "22", 0.0, 1000, "AC", "A"),  # deletion
        ("rs_indel2", "22", 0.0, 2000, "A", "AT"),  # insertion
        ("rs_snp", "22", 0.0, 3000, "A", "G"),
    ])
    counters = build_sites_vcf(df, tmp_path / "out.vcf", aadr_build="hg19")
    assert counters.non_snp_drops == 2
    assert counters.rows_written == 1


def test_non_snp_filter_drops_non_acgt(tmp_path: Path) -> None:
    """Non-ACGT alleles (N, ambiguity codes) drop."""
    df = _df_from_rows([
        ("rs_n", "22", 0.0, 1000, "A", "N"),
        ("rs_w", "22", 0.0, 2000, "W", "S"),
    ])
    counters = build_sites_vcf(df, tmp_path / "out.vcf", aadr_build="hg19")
    assert counters.non_snp_drops == 2
    assert counters.rows_written == 0


# --- Chrom normalization + sort order ---


def test_unrecognized_chroms_silently_skipped(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Alt-contigs / decoys / invalid chrom encodings are skipped with DEBUG logging."""
    import logging
    caplog.set_level(logging.DEBUG, logger="pileup_aadr.sites_vcf")
    df = _df_from_rows([
        ("rs_alt", "chr1_KI270706v1_random", 0.0, 1000, "A", "G"),
        ("rs_unknown", "999", 0.0, 2000, "C", "T"),
        ("rs_ok", "22", 0.0, 3000, "A", "G"),
    ])
    counters = build_sites_vcf(df, tmp_path / "out.vcf", aadr_build="hg19")
    assert counters.rows_written == 1
    # Counters don't increment for unrecognized chroms (they're not "non-SNP" and not
    # "palindromic" — they're just skipped). Verified via the rows_written delta.
    assert counters.non_snp_drops == 0
    assert counters.palindrome_drops == 0


def test_rows_preserve_chrom_order_from_parse(tmp_path: Path) -> None:
    """build_sites_vcf preserves input row order (sort contract is in parse_aadr_snp).

    Input is pre-sorted as parse_aadr_snp guarantees (CHROM_ORDER x pos_bp); output
    must preserve that order. chr2 before chr22 before chrX per CHROM_ORDER.
    """
    df = _df_from_rows([
        # pre-sorted: chr2 < chr22 < chrX, and within chr22 early < late
        ("rs_2_only", "2", 0.0, 5000, "A", "G"),
        ("rs_22_early", "22", 0.0, 1000, "A", "G"),
        ("rs_22_late", "22", 0.0, 9000, "A", "G"),
        ("rs_x_only", "23", 0.0, 7000, "A", "G"),
    ])
    out_vcf = tmp_path / "out.vcf"
    build_sites_vcf(df, out_vcf, aadr_build="hg19")

    data_lines = [
        line.split("\t")[2]  # ID column
        for line in out_vcf.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert data_lines == ["rs_2_only", "rs_22_early", "rs_22_late", "rs_x_only"]


# --- aadr_build switch ---


def test_aadr_build_hg38_switches_chrom_lengths_table(tmp_path: Path) -> None:
    """When aadr_build='hg38', ##contig lengths come from HG38_CHROM_LENGTHS."""
    from pileup_aadr.chrom_lengths import HG38_CHROM_LENGTHS

    df = _df_from_rows([
        ("rs1", "22", 0.0, 1000, "A", "G"),
    ])
    out_vcf = tmp_path / "out.vcf"
    build_sites_vcf(df, out_vcf, aadr_build="hg38")
    content = out_vcf.read_text()
    expected = HG38_CHROM_LENGTHS["chr22"]
    assert f"##contig=<ID=chr22,length={expected}>" in content


def test_palindrome_count_excludes_indel_palindromes(tmp_path: Path) -> None:
    """Indel-palindromes (AT/TA as len>1 or non-ACGT) count as non_snp_drops, not palindrome_drops."""
    df = _df_from_rows([
        # Genuine indel that looks palindromic — should be non_snp_drop, not palindrome_drop
        ("rs_at_indel", "22", 0.0, 1000, "AT", "TA"),
        # Genuine SNP palindrome — should be palindrome_drop
        ("rs_at_snp", "22", 0.0, 2000, "A", "T"),
        # Clean SNP — should pass through
        ("rs_ag_snp", "22", 0.0, 3000, "A", "G"),
    ])
    counters = build_sites_vcf(df, tmp_path / "out.vcf", aadr_build="hg19")
    assert counters.non_snp_drops == 1   # only the indel
    assert counters.palindrome_drops == 1  # only the genuine SNP palindrome
    assert counters.rows_written == 1


def test_vectorized_filter_and_sort_multi_chrom(tmp_path: Path) -> None:
    """Multi-chrom mix of indels/palindromes/valid rows: correct counts and CHROM_ORDER output.

    Input is pre-sorted as parse_aadr_snp guarantees (CHROM_ORDER x pos_bp within chrom).
    """
    df = _df_from_rows([
        # pre-sorted: chr1 < chr2 < chr22
        ("rs_1_indel", "1", 0.0, 200, "AC", "A"),   # non-SNP
        ("rs_1_snp", "1", 0.0, 300, "C", "T"),      # valid
        ("rs_2_snp", "2", 0.0, 400, "G", "A"),      # valid
        ("rs_22_snp", "22", 0.0, 100, "A", "G"),    # valid (early pos)
        ("rs_22_pal", "22", 0.0, 500, "C", "G"),    # palindrome
    ])
    out_vcf = tmp_path / "out.vcf"
    counters = build_sites_vcf(df, out_vcf, aadr_build="hg19")

    assert counters.non_snp_drops == 1
    assert counters.palindrome_drops == 1
    assert counters.rows_written == 3

    ids = [
        line.split("\t")[2]
        for line in out_vcf.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    # chr1 < chr2 < chr22 in CHROM_ORDER; within chr22 pos 100 only
    assert ids == ["rs_1_snp", "rs_2_snp", "rs_22_snp"]


def test_empty_input_writes_header_only(tmp_path: Path) -> None:
    """Empty AADR DataFrame → VCF with header but no data rows."""
    df = _df_from_rows([])
    out_vcf = tmp_path / "out.vcf"
    counters = build_sites_vcf(df, out_vcf, aadr_build="hg19")
    assert counters.rows_written == 0
    content = out_vcf.read_text()
    assert "##fileformat=VCFv4.2" in content
    assert "#CHROM\tPOS" in content
    # No data rows
    data_lines = [
        line for line in content.splitlines()
        if line and not line.startswith("#")
    ]
    assert data_lines == []
