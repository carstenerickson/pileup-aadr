"""Tests for format_detect.py — AADR .snp parser, build detection, normalize_chrom."""
from __future__ import annotations

from pathlib import Path

import pytest

from pileup_aadr.errors import (
    AADRDuplicateRsidError,
    AADRParseError,
    UnsupportedAADRBuild,
    UnsupportedReferenceBuild,
)
from pileup_aadr.format_detect import (
    detect_aadr_build,
    detect_bam_build,
    normalize_chrom,
    parse_aadr_snp,
)

# --- normalize_chrom ---


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("1", "chr1"),
        ("chr1", "chr1"),
        ("22", "chr22"),
        ("X", "chrX"),
        ("23", "chrX"),
        ("chrX", "chrX"),
        ("Y", "chrY"),
        ("24", "chrY"),
        ("MT", "chrM"),
        ("chrMT", "chrM"),
        ("M", "chrM"),
        ("90", "chrM"),
        ("91", "chrXY"),
    ],
)
def test_normalize_chrom_known_values(raw: str, expected: str) -> None:
    assert normalize_chrom(raw) == expected


@pytest.mark.parametrize(
    "raw",
    ["chr1_KI270706v1_random", "chrUn_GL000218v1", "chrEBV", "decoy", "999", ""],
)
def test_normalize_chrom_returns_none_for_alt_contigs(raw: str) -> None:
    assert normalize_chrom(raw) is None


# --- parse_aadr_snp ---


def test_parse_aadr_snp_chr22_slice(aadr_chr22_slice: Path) -> None:
    """The committed slice parses cleanly and indexes by rsid."""
    df = parse_aadr_snp(aadr_chr22_slice)
    assert len(df) == 50
    assert df.index.name == "rsid"
    assert df.index.is_unique
    # Spot-check known row
    assert df.loc["rs17428495", "chrom_int"] == "22"
    assert df.loc["rs17428495", "ref"] == "G"
    assert df.loc["rs17428495", "alt"] == "A"
    assert df.loc["rs17428495", "pos_bp"] == 16071624


def test_parse_aadr_snp_handles_multi_space_padding(tmp_path: Path) -> None:
    """AADR's whitespace is multi-space-padded for column alignment, NOT tab-separated.
    str.split() handles both gracefully."""
    snp = tmp_path / "padded.snp"
    snp.write_text(
        # Mix of multi-space, tab, single-space — all valid per AADR
        "          rs1     22     0.0     1000  A G\n"
        "rs2\t22\t0.0\t2000\tC\tT\n"
        "rs3 22 0.0 3000 G C\n"
    )
    df = parse_aadr_snp(snp)
    assert len(df) == 3
    assert list(df.index) == ["rs1", "rs2", "rs3"]


def test_parse_aadr_snp_skips_comments_and_blanks(tmp_path: Path) -> None:
    snp = tmp_path / "with_comments.snp"
    snp.write_text(
        "# This is a comment\n"
        "\n"
        "rs1     22     0.0     1000  A G\n"
        "# Another comment\n"
        "rs2     22     0.0     2000  C T\n"
    )
    df = parse_aadr_snp(snp)
    assert len(df) == 2


def test_parse_aadr_snp_rejects_duplicate_rsid(tmp_path: Path) -> None:
    snp = tmp_path / "dup.snp"
    snp.write_text(
        "rs1     22     0.0     1000  A G\n"
        "rs2     22     0.0     2000  C T\n"
        "rs1     22     0.0     3000  G A\n"  # duplicate
    )
    with pytest.raises(AADRDuplicateRsidError, match="rs1"):
        parse_aadr_snp(snp)


def test_parse_aadr_snp_rejects_wrong_column_count(tmp_path: Path) -> None:
    snp = tmp_path / "wrong_cols.snp"
    snp.write_text("rs1 22 0.0 1000 A\n")  # only 5 columns
    with pytest.raises(AADRParseError, match="6 columns"):
        parse_aadr_snp(snp)


def test_parse_aadr_snp_rejects_non_acgt_alleles(tmp_path: Path) -> None:
    snp = tmp_path / "bad_allele.snp"
    snp.write_text("rs1 22 0.0 1000 A N\n")  # N is not ACGT
    with pytest.raises(AADRParseError, match="non-ACGT"):
        parse_aadr_snp(snp)


def test_parse_aadr_snp_rejects_unparseable_position(tmp_path: Path) -> None:
    snp = tmp_path / "bad_pos.snp"
    snp.write_text("rs1 22 0.0 NOT_A_NUMBER A G\n")
    with pytest.raises(AADRParseError):
        parse_aadr_snp(snp)


# --- detect_aadr_build ---


def test_detect_aadr_build_override_short_circuits(aadr_chr22_slice: Path) -> None:
    """Explicit override returns immediately without consulting data."""
    df = parse_aadr_snp(aadr_chr22_slice)
    assert detect_aadr_build(df, override="hg19") == "hg19"
    assert detect_aadr_build(df, override="hg38") == "hg38"


def test_detect_aadr_build_no_chr1_rows_raises(aadr_chr22_slice: Path) -> None:
    """The chr22 slice has no chr1 rows — auto-detection has nothing to consult."""
    df = parse_aadr_snp(aadr_chr22_slice)
    with pytest.raises(UnsupportedAADRBuild, match="chr1"):
        detect_aadr_build(df, override="auto")


def test_detect_aadr_build_hg19(tmp_path: Path) -> None:
    """A row with chr1 max position closer to hg19 chr1 length than hg38's → hg19.

    HG19_CHR1_LENGTH = 249,250,621
    HG38_CHR1_LENGTH = 248,956,422 (294 KB shorter)
    Pick a position ≤ 100K from hg19's end so it's unambiguously closer to hg19.
    """
    snp = tmp_path / "hg19_chr1.snp"
    snp.write_text("rs1 1 0.0 249200000 A G\n")  # 50,621 from hg19; 243,578 from hg38
    df = parse_aadr_snp(snp)
    assert detect_aadr_build(df, override="auto") == "hg19"


def test_detect_aadr_build_hg38(tmp_path: Path) -> None:
    """A row with chr1 max position closer to hg38 chr1 length than hg19's → hg38."""
    snp = tmp_path / "hg38_chr1.snp"
    snp.write_text("rs1 1 0.0 248900000 A G\n")  # 56,422 from hg38; 350,621 from hg19
    df = parse_aadr_snp(snp)
    assert detect_aadr_build(df, override="auto") == "hg38"


# --- detect_bam_build (issue #1 regression) ---


def _make_bam_with_chr1_length(path: Path, chr1_length: int) -> Path:
    """Write a header-only BAM with a configurable chr1 @SQ length."""
    import pysam
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [
            {"SN": "chr1", "LN": chr1_length},
            {"SN": "chr22", "LN": 51_304_566},
        ],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass
    return path


def test_detect_bam_build_hg38_not_misclassified_as_hg19(tmp_path: Path) -> None:
    """Issue #1 regression: hg38 chr1 (248,956,422) is only 294 KB shorter than
    hg19 (249,250,621) — well INSIDE the ±1 Mb tolerance window. The pre-fix
    first-match-wins logic always returned 'hg19' for hg38 BAMs because the
    hg19 check fired first within tolerance. Closest-match must pick hg38.
    """
    bam = _make_bam_with_chr1_length(tmp_path / "hg38.bam", 248_956_422)
    assert detect_bam_build(bam, override="auto") == "hg38"


def test_detect_bam_build_hg19(tmp_path: Path) -> None:
    """Sanity: hg19 chr1 (249,250,621) detects as hg19."""
    bam = _make_bam_with_chr1_length(tmp_path / "hg19.bam", 249_250_621)
    assert detect_bam_build(bam, override="auto") == "hg19"


def test_detect_bam_build_override_short_circuits(tmp_path: Path) -> None:
    """Explicit override returns immediately without consulting @SQ."""
    bam = _make_bam_with_chr1_length(tmp_path / "any.bam", 248_956_422)
    assert detect_bam_build(bam, override="hg19") == "hg19"
    assert detect_bam_build(bam, override="hg38") == "hg38"


def test_detect_bam_build_unknown_assembly_raises(tmp_path: Path) -> None:
    """T2T-CHM13 chr1 (248_387_328) is >500 KB from hg38 — outside tolerance."""
    bam = _make_bam_with_chr1_length(tmp_path / "t2t.bam", 200_000_000)
    with pytest.raises(UnsupportedReferenceBuild, match=r"neither hg19 .*nor hg38"):
        detect_bam_build(bam, override="auto")
