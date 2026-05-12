"""Tests for ref_resolve.py — BAM @PG FASTA extraction + build verification."""
from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from pileup_aadr.errors import (
    ReferenceFastaBuildMismatchError,
    ReferenceFastaNotFound,
)
from pileup_aadr.format_detect import HG19_CHR1_LENGTH, HG38_CHR1_LENGTH
from pileup_aadr.ref_resolve import (
    _extract_ref_from_bam_pg,
    resolve_ref_fasta,
    verify_fasta_matches_bam_build,
)


def _write_fai(fasta_path: Path, chr1_length: int) -> None:
    """Write a minimal .fai sidecar with chr1 having the given length."""
    fai = Path(f"{fasta_path}.fai")
    fai.write_text(f"chr1\t{chr1_length}\t6\t60\t61\n")


# --- verify_fasta_matches_bam_build ---


def test_verify_fasta_matches_bam_build_hg38_ok(tmp_path: Path) -> None:
    """hg38 BAM + hg38 FASTA (chr1 length within ±1 Mb) → no raise."""
    fasta = tmp_path / "hg38.fa"
    fasta.touch()
    _write_fai(fasta, HG38_CHR1_LENGTH)
    verify_fasta_matches_bam_build(fasta, "hg38")  # no raise


def test_verify_fasta_matches_bam_build_hg19_ok(tmp_path: Path) -> None:
    fasta = tmp_path / "hg19.fa"
    fasta.touch()
    _write_fai(fasta, HG19_CHR1_LENGTH)
    verify_fasta_matches_bam_build(fasta, "hg19")  # no raise


def test_verify_fasta_matches_bam_build_mismatch_raises(tmp_path: Path) -> None:
    """hg38 BAM + hg19 FASTA → ReferenceFastaBuildMismatchError."""
    fasta = tmp_path / "wrong.fa"
    fasta.touch()
    _write_fai(fasta, HG19_CHR1_LENGTH)  # FASTA is hg19...
    with pytest.raises(ReferenceFastaBuildMismatchError, match=r"hg19.*hg38"):
        verify_fasta_matches_bam_build(fasta, "hg38")  # ...but BAM is hg38


def test_verify_fasta_matches_bam_build_skips_when_fai_missing(tmp_path: Path) -> None:
    """No .fai sidecar → skip verification (Picard will discover the issue itself)."""
    fasta = tmp_path / "no_fai.fa"
    fasta.touch()
    # No .fai created
    verify_fasta_matches_bam_build(fasta, "hg38")  # no raise; just logs at DEBUG


def test_verify_fasta_matches_bam_build_skips_when_no_chr1_in_fai(tmp_path: Path) -> None:
    """`.fai` exists but lists no chr1 — skip (defensive; pysam tolerates oddly-indexed FASTAs)."""
    fasta = tmp_path / "no_chr1.fa"
    fasta.touch()
    fai = Path(f"{fasta}.fai")
    fai.write_text("chr2\t100000\t6\t60\t61\nchr3\t200000\t6\t60\t61\n")
    verify_fasta_matches_bam_build(fasta, "hg38")  # no raise


def test_verify_fasta_matches_bam_build_within_tolerance(tmp_path: Path) -> None:
    """Small patch-level drift accepted as long as it stays closest to the claimed build.

    hg19 and hg38 chr1 differ by only 294 KB, so the tolerance window for "still
    closest to hg38" is narrow. 100 KB drift toward hg19 keeps hg38 closest
    (294-100=194 KB > 100 KB). Larger drifts cross into the hg19 range.
    """
    fasta = tmp_path / "patched.fa"
    fasta.touch()
    _write_fai(fasta, HG38_CHR1_LENGTH - 100_000)  # 100 KB drift away from hg19
    verify_fasta_matches_bam_build(fasta, "hg38")  # no raise — still closest to hg38


# --- _extract_ref_from_bam_pg ---


@pytest.fixture
def bam_with_pg_path(tmp_path: Path) -> Path:
    """Construct a BAM whose @PG line embeds a FASTA path that exists on disk."""
    fasta = tmp_path / "ref.fa"
    fasta.touch()
    _write_fai(fasta, HG38_CHR1_LENGTH)
    bam = tmp_path / "with_pg.bam"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": HG38_CHR1_LENGTH}],
        "PG": [{"ID": "bwa", "PN": "bwa", "CL": f"bwa mem -R '@RG\\tID:test' {fasta} fq1"}],
    }
    with pysam.AlignmentFile(str(bam), "wb", header=header):
        pass
    pysam.sort("-o", str(bam), str(bam))
    pysam.index(str(bam))
    return bam


def test_extract_ref_from_bam_pg_finds_existing_fasta(bam_with_pg_path: Path) -> None:
    """The @PG line's FASTA path exists on disk → return it."""
    extracted = _extract_ref_from_bam_pg(bam_with_pg_path)
    assert extracted.exists()
    assert extracted.suffix == ".fa"


def test_extract_ref_from_bam_pg_raises_when_no_fasta_in_pg(tmp_path: Path) -> None:
    """BAM with no @PG records → ReferenceFastaNotFound."""
    bam = tmp_path / "no_pg.bam"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": HG38_CHR1_LENGTH}],
    }
    with pysam.AlignmentFile(str(bam), "wb", header=header):
        pass
    pysam.sort("-o", str(bam), str(bam))
    pysam.index(str(bam))
    with pytest.raises(ReferenceFastaNotFound, match=r"no readable \.fa path"):
        _extract_ref_from_bam_pg(bam)


def test_extract_ref_from_bam_pg_raises_when_pg_path_missing(tmp_path: Path) -> None:
    """@PG embeds a path but the file doesn't exist → raise."""
    bam = tmp_path / "stale_pg.bam"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": HG38_CHR1_LENGTH}],
        "PG": [{"ID": "bwa", "CL": "bwa mem /nonexistent/path/to/ref.fa fq1"}],
    }
    with pysam.AlignmentFile(str(bam), "wb", header=header):
        pass
    pysam.sort("-o", str(bam), str(bam))
    pysam.index(str(bam))
    with pytest.raises(ReferenceFastaNotFound):
        _extract_ref_from_bam_pg(bam)


# --- resolve_ref_fasta (full 3-tier) ---


def test_resolve_ref_fasta_uses_cli_explicit(tmp_path: Path, bam_with_pg_path: Path) -> None:
    """Explicit --ref-fasta wins, but is still build-verified."""
    fasta = tmp_path / "explicit.fa"
    fasta.touch()
    _write_fai(fasta, HG38_CHR1_LENGTH)
    result = resolve_ref_fasta(cli_ref=fasta, bam=bam_with_pg_path, bam_build="hg38")
    assert result == fasta


def test_resolve_ref_fasta_explicit_raises_on_build_mismatch(
    tmp_path: Path, bam_with_pg_path: Path
) -> None:
    """Explicit --ref-fasta is NOT exempt from build verification (LLD H11)."""
    fasta = tmp_path / "wrong_build.fa"
    fasta.touch()
    _write_fai(fasta, HG19_CHR1_LENGTH)
    with pytest.raises(ReferenceFastaBuildMismatchError):
        resolve_ref_fasta(cli_ref=fasta, bam=bam_with_pg_path, bam_build="hg38")


def test_resolve_ref_fasta_uses_env_dir(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, bam_with_pg_path: Path
) -> None:
    """$PILEUP_AADR_REF_DIR/<build>.fa picked up when --ref-fasta not given."""
    env_dir = tmp_path / "refs"
    env_dir.mkdir()
    fasta = env_dir / "hg38.fa"
    fasta.touch()
    _write_fai(fasta, HG38_CHR1_LENGTH)
    monkeypatch.setenv("PILEUP_AADR_REF_DIR", str(env_dir))
    result = resolve_ref_fasta(cli_ref=None, bam=bam_with_pg_path, bam_build="hg38")
    assert result == fasta


def test_resolve_ref_fasta_falls_through_to_bam_pg(
    monkeypatch: pytest.MonkeyPatch, bam_with_pg_path: Path
) -> None:
    """No --ref-fasta, no env → @PG extraction returns the BAM-recorded path."""
    monkeypatch.delenv("PILEUP_AADR_REF_DIR", raising=False)
    result = resolve_ref_fasta(cli_ref=None, bam=bam_with_pg_path, bam_build="hg38")
    assert result.exists()
    assert result.suffix == ".fa"
