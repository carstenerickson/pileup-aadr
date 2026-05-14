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


# --- chr20 fallback (v0.2 enhancement; customer-suggested) ---


def _make_bam_with_chr20_only(path: Path, chr20_length: int) -> Path:
    """A BAM @SQ with chr20 but no chr1 — exercises the fallback path."""
    import pysam
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr20", "LN": chr20_length}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass
    return path


def test_detect_bam_build_chr20_fallback_hg19(tmp_path: Path) -> None:
    """When chr1 is absent, fall back to chr20. hg19 chr20 = 63,025,520."""
    bam = _make_bam_with_chr20_only(tmp_path / "hg19_chr20.bam", 63_025_520)
    assert detect_bam_build(bam, override="auto") == "hg19"


def test_detect_bam_build_chr20_fallback_hg38(tmp_path: Path) -> None:
    """When chr1 is absent, fall back to chr20. hg38 chr20 = 64,444,167."""
    bam = _make_bam_with_chr20_only(tmp_path / "hg38_chr20.bam", 64_444_167)
    assert detect_bam_build(bam, override="auto") == "hg38"


def test_detect_bam_build_no_chr1_no_chr20_raises(tmp_path: Path) -> None:
    """BAM with only chrY (or any non-chr1/chr20 chrom) — raises with diagnostic."""
    import pysam
    bam = tmp_path / "chrY_only.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrY", "LN": 57_227_415}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header):
        pass
    with pytest.raises(UnsupportedReferenceBuild, match=r"chr1.*nor chr20"):
        detect_bam_build(bam, override="auto")


def test_detect_aadr_build_chr20_fallback_hg19(tmp_path: Path) -> None:
    """When AADR has no chr1 rows, fall back to chr20. hg19 chr20 ≈ 63.0 Mb."""
    snp = tmp_path / "hg19_chr20.snp"
    snp.write_text("rs1 20 0.0 62900000 A G\n")  # within 5 Mb of hg19 chr20 end
    df = parse_aadr_snp(snp)
    assert detect_aadr_build(df, override="auto") == "hg19"


def test_detect_aadr_build_chr20_fallback_hg38(tmp_path: Path) -> None:
    """When AADR has no chr1 rows, fall back to chr20. hg38 chr20 ≈ 64.4 Mb."""
    snp = tmp_path / "hg38_chr20.snp"
    snp.write_text("rs1 20 0.0 64400000 A G\n")  # within 5 Mb of hg38 chr20 end
    df = parse_aadr_snp(snp)
    assert detect_aadr_build(df, override="auto") == "hg38"


def test_detect_aadr_build_no_chr1_no_chr20_raises(tmp_path: Path) -> None:
    """AADR slice with only chr22 (e.g., the existing chr22 fixture) — diagnostic
    names BOTH anchors (chr1 + chr20) so users can see what was tried."""
    snp = tmp_path / "chr22_only.snp"
    snp.write_text("rs1 22 0.0 51000000 A G\n")
    df = parse_aadr_snp(snp)
    with pytest.raises(UnsupportedAADRBuild, match=r"chr1.*OR chr20"):
        detect_aadr_build(df, override="auto")


# --- classify_aadr_chrom_set (v0.2 panel-class enhancement) ---


def _aadr_df_from_chroms(tmp_path: Path, chrom_ints: list[str]) -> pd.DataFrame:  # noqa: F821
    """Build a tiny AADR DataFrame with one row per requested chrom_int."""
    snp = tmp_path / "tmp.snp"
    snp.write_text("\n".join(
        f"rs{i} {c} 0.0 {1000 * (i + 1)} A G" for i, c in enumerate(chrom_ints)
    ) + "\n")
    return parse_aadr_snp(snp)


def test_classify_full_panel_is_autosomes_plus_sex(tmp_path: Path) -> None:
    """A panel with chr1-22 + chr X + Y (typical 1240k / HO shape)."""
    from pileup_aadr.format_detect import classify_aadr_chrom_set
    chroms = [str(i) for i in range(1, 23)] + ["23", "24"]
    df = _aadr_df_from_chroms(tmp_path, chroms)
    assert classify_aadr_chrom_set(df) == "autosomes+sex"


def test_classify_autosomes_only(tmp_path: Path) -> None:
    """A 1240k autosomes-only subset: all 22 autosomes, no sex/MT."""
    from pileup_aadr.format_detect import classify_aadr_chrom_set
    df = _aadr_df_from_chroms(tmp_path, [str(i) for i in range(1, 23)])
    assert classify_aadr_chrom_set(df) == "autosomes_only"


def test_classify_chry_only(tmp_path: Path) -> None:
    """chrY-only AADR slice (haplogroup workflows)."""
    from pileup_aadr.format_detect import classify_aadr_chrom_set
    df = _aadr_df_from_chroms(tmp_path, ["24"])
    assert classify_aadr_chrom_set(df) == "chrY_only"


def test_classify_chrm_only(tmp_path: Path) -> None:
    """chrM-only AADR slice (mtDNA workflows)."""
    from pileup_aadr.format_detect import classify_aadr_chrom_set
    df = _aadr_df_from_chroms(tmp_path, ["90"])
    assert classify_aadr_chrom_set(df) == "chrM_only"


def test_classify_sex_only(tmp_path: Path) -> None:
    """chrX + chrY only (sex-chromosome-specific panel)."""
    from pileup_aadr.format_detect import classify_aadr_chrom_set
    df = _aadr_df_from_chroms(tmp_path, ["23", "24"])
    assert classify_aadr_chrom_set(df) == "sex_only"


def test_classify_chr22_only_is_custom(tmp_path: Path) -> None:
    """Single-autosome slice (the existing chr22 test fixture shape) is `custom`."""
    from pileup_aadr.format_detect import classify_aadr_chrom_set
    df = _aadr_df_from_chroms(tmp_path, ["22"])
    assert classify_aadr_chrom_set(df) == "custom"


# --- C1-C4: AADR .snp parse cache ---

def _write_minimal_snp(path: Path) -> None:
    """Write a minimal valid AADR .snp with 2 rows."""
    path.write_text(
        "rs1\t1\t0.0\t1000\tA\tG\n"
        "rs2\t1\t0.0\t2000\tC\tT\n"
    )


def test_cache_miss_populates_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """C1: first parse writes feather; result equals second parse (cache hit)."""
    snp = tmp_path / "sites.snp"
    _write_minimal_snp(snp)
    cache_dir = tmp_path / "cache"
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache_dir))
    monkeypatch.delenv("PILEUP_AADR_DISABLE_SNP_CACHE", raising=False)

    df1 = parse_aadr_snp(snp)
    feather_files = list((cache_dir / "pileup-aadr" / "snp").glob("*.feather"))
    assert len(feather_files) == 1, "cache write should produce exactly one feather"

    df2 = parse_aadr_snp(snp)
    assert list(df1.index) == list(df2.index)
    assert list(df1["pos_bp"]) == list(df2["pos_bp"])


def test_cache_hit_skips_full_parse(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """C2: second call uses pd.read_feather (cache hit path), not the text parser."""
    import unittest.mock as mock
    import pandas as pd
    import pileup_aadr.format_detect as fd

    snp = tmp_path / "sites.snp"
    _write_minimal_snp(snp)
    cache_dir = tmp_path / "cache"
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache_dir))
    monkeypatch.delenv("PILEUP_AADR_DISABLE_SNP_CACHE", raising=False)

    parse_aadr_snp(snp)  # cold parse → writes cache

    with mock.patch.object(fd.pd, "read_feather", wraps=pd.read_feather) as spy:
        df2 = parse_aadr_snp(snp)
        assert spy.call_count == 1, "cache hit must call read_feather once"

    assert "rs1" in df2.index


def test_cache_invalidated_on_content_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """C3: different file content → different SHA256 → cache miss → new entry."""
    cache_dir = tmp_path / "cache"
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache_dir))
    monkeypatch.delenv("PILEUP_AADR_DISABLE_SNP_CACHE", raising=False)

    snp = tmp_path / "sites.snp"
    _write_minimal_snp(snp)
    parse_aadr_snp(snp)

    snp.write_text(
        "rs1\t1\t0.0\t1000\tA\tG\n"
        "rs2\t1\t0.0\t2000\tC\tT\n"
        "rs3\t2\t0.0\t3000\tA\tC\n"
    )
    df2 = parse_aadr_snp(snp)
    assert "rs3" in df2.index

    feather_files = list((cache_dir / "pileup-aadr" / "snp").glob("*.feather"))
    assert len(feather_files) == 2, "two distinct content hashes → two cache entries"


def test_cache_disabled_env_var_skips_read_and_write(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """C4: PILEUP_AADR_DISABLE_SNP_CACHE=1 → no feather written; re-parse always re-reads."""
    cache_dir = tmp_path / "cache"
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache_dir))
    monkeypatch.setenv("PILEUP_AADR_DISABLE_SNP_CACHE", "1")

    snp = tmp_path / "sites.snp"
    _write_minimal_snp(snp)
    parse_aadr_snp(snp)

    feather_files = list((cache_dir / "pileup-aadr" / "snp").glob("*.feather")) if (
        cache_dir / "pileup-aadr" / "snp"
    ).exists() else []
    assert feather_files == [], "cache write must be skipped when env var is set"
