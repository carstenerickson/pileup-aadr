"""Tests for Stage 1 — parse_picard_stderr + parse_rejected_vcf + lift_aadr_sites.

Most tests use captured Picard stderr fixtures from v2.1 verification (real Picard
3.3.0 output on `ancestrytracke-f`). The lift_aadr_sites integration test mocks
ToolWrapper.run to avoid needing Picard in the test env.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from pileup_aadr.counters import Stage1InputFilters
from pileup_aadr.errors import LiftoverYieldError, PileupAadrInternalError
from pileup_aadr.lift import (
    lift_aadr_sites,
    parse_picard_stderr,
    parse_rejected_vcf,
)

FIXTURES_STDERR = Path(__file__).parent / "fixtures" / "stderr"


# --- parse_picard_stderr ---


def test_parse_picard_stderr_clean_run() -> None:
    """The 100-site clean run from v2.1: 100 processed, 0 failed, 0 swapped."""
    text = (FIXTURES_STDERR / "picard_clean.stderr").read_text()
    parsed = parse_picard_stderr(text)
    assert parsed["processed_count"] == 100
    assert parsed["failed_count"] == 0
    assert parsed["mismatched_ref_count"] == 0
    assert parsed["lifted_count"] == 100
    assert parsed["swapped_count"] == 0


def test_parse_picard_stderr_partial_yield() -> None:
    """The 5000-site v2.1 run: 5000 processed, 2 failed, 16 swapped."""
    text = (FIXTURES_STDERR / "picard_partial_yield.stderr").read_text()
    parsed = parse_picard_stderr(text)
    assert parsed["processed_count"] == 5000
    assert parsed["failed_count"] == 2
    assert parsed["lifted_count"] == 4998
    assert parsed["swapped_count"] == 16


def test_parse_picard_stderr_swap_line_optional() -> None:
    """If Picard ever omits the 'X variants were lifted by swapping' line, default to 0."""
    text = """
    INFO	LiftoverVcf	Processed 100 variants.
    INFO	LiftoverVcf	0 variants failed to liftover.
    INFO	LiftoverVcf	0 variants lifted over but had mismatching reference alleles after lift over.
    INFO	LiftoverVcf	0.0000% of variants were not successfully lifted over
    """
    parsed = parse_picard_stderr(text)
    assert parsed["swapped_count"] == 0


def test_parse_picard_stderr_missing_required_pattern_raises() -> None:
    """Missing a REQUIRED stderr line (e.g., 'Processed N variants') → PileupAadrInternalError."""
    text = "Some unrelated output that doesn't match Picard's structured stderr."
    with pytest.raises(PileupAadrInternalError, match="processed"):
        parse_picard_stderr(text)


# --- parse_rejected_vcf ---


def _write_rejected_vcf(path: Path, records: list[tuple[str, int, str]]) -> None:
    """Build a minimal rejected.vcf with (chrom, pos, FILTER) tuples."""
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=249250621>",
        '##FILTER=<ID=NoTarget,Description="No chain alignment">',
        '##FILTER=<ID=MismatchedRefAllele,Description="REF allele mismatch">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for chrom, pos, filt in records:
        lines.append(f"{chrom}\t{pos}\trs1\tA\tG\t.\t{filt}\t.")
    path.write_text("\n".join(lines) + "\n")


def test_parse_rejected_vcf_empty(tmp_path: Path) -> None:
    """No rejected.vcf at the path → all-zero counts (defensive)."""
    counts = parse_rejected_vcf(tmp_path / "missing.vcf")
    assert counts == {
        "NoTarget": 0,
        "MismatchedRefAllele": 0,
        "IndelStraddlesMultipleIntervals": 0,
        "SwappedAlleles": 0,
        "other": 0,
    }


def test_parse_rejected_vcf_categorizes_filters(tmp_path: Path) -> None:
    """Each known FILTER value increments its bucket; unknowns go to 'other'."""
    rej = tmp_path / "rejected.vcf"
    _write_rejected_vcf(
        rej,
        [
            ("chr1", 1000, "NoTarget"),
            ("chr1", 2000, "NoTarget"),
            ("chr1", 3000, "MismatchedRefAllele"),
            ("chr1", 4000, "WeirdUnknownFilter"),
        ],
    )
    counts = parse_rejected_vcf(rej)
    assert counts["NoTarget"] == 2
    assert counts["MismatchedRefAllele"] == 1
    assert counts["other"] == 1
    assert counts["IndelStraddlesMultipleIntervals"] == 0


# --- lift_aadr_sites integration (mocked Picard subprocess) ---


@pytest.fixture
def mock_picard_lift(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """Mock ToolWrapper.run to write a fake Picard stderr + empty rejected.vcf,
    leaving the OUTPUT VCF up to the caller to either pre-create or ignore.
    Returns the tmp_path for callers to set up files in.
    """
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))

    # Mock subprocess.run for the version probe (called from ToolWrapper.__init__)
    import subprocess
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    # Mock ToolWrapper.run to write the captured-stderr fixture + return success
    from pileup_aadr import lift, tool_wrapper

    def fake_run(self: object, args: list[str], *, capture_stderr_to: Path, **_kw: object) -> object:
        # Copy the v2.1 captured Picard stderr into the path the caller expects
        capture_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        capture_stderr_to.write_text(
            (FIXTURES_STDERR / "picard_clean.stderr").read_text()
        )
        return tool_wrapper.ToolRunResult(
            exit_code=0,
            stdout=None,
            stderr_path=capture_stderr_to,
            stderr_text=None,
            wallclock_seconds=0.1,
            peak_rss_mb=None,
        )

    monkeypatch.setattr(lift.ToolWrapper, "run", fake_run)
    monkeypatch.setattr(lift.ToolWrapper, "_check_version", lambda _self: None)
    return tmp_path


def test_lift_aadr_sites_clean_run(mock_picard_lift: Path) -> None:
    """100/100 sites lifted (per the captured stderr) → yield = 100% → PASS gate."""
    sites_vcf = mock_picard_lift / "sites.vcf"
    sites_vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    fasta = mock_picard_lift / "ref.fa"
    fasta.touch()
    chain = mock_picard_lift / "chain.gz"
    chain.touch()

    counters = lift_aadr_sites(
        sites_vcf_path=sites_vcf,
        chain_path=chain,
        target_fasta_path=fasta,
        output_lifted_vcf=mock_picard_lift / "lifted.vcf",
        output_rejected_vcf=mock_picard_lift / "rejected.vcf",
        input_filter_counters=Stage1InputFilters(
            palindrome_drops=0, non_snp_drops=0, non_autosome_drops=0, rows_written=100
        ),
    )

    assert counters.input_sites_after_filters == 100
    assert counters.lifted_sites == 100
    assert counters.liftover_yield_pct == 100.0
    assert counters.liftover_yield_warning is False
    assert counters.swapped_alleles_count == 0


def test_lift_aadr_sites_low_yield_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Synthetic stderr showing 50% yield → LiftoverYieldError below default 70% gate."""
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))

    import subprocess
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    fake_low_yield_stderr = """
    INFO	LiftoverVcf	Processed 100 variants.
    INFO	LiftoverVcf	50 variants failed to liftover.
    INFO	LiftoverVcf	0 variants lifted over but had mismatching reference alleles after lift over.
    INFO	LiftoverVcf	50.0000% of variants were not successfully lifted over
    INFO	LiftoverVcf	0 variants were lifted by swapping REF/ALT alleles.
    """

    from pileup_aadr import lift, tool_wrapper

    def fake_run(self: object, args: list[str], *, capture_stderr_to: Path, **_kw: object) -> object:
        capture_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        capture_stderr_to.write_text(fake_low_yield_stderr)
        return tool_wrapper.ToolRunResult(
            exit_code=0, stdout=None, stderr_path=capture_stderr_to,
            stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None,
        )

    monkeypatch.setattr(lift.ToolWrapper, "run", fake_run)
    monkeypatch.setattr(lift.ToolWrapper, "_check_version", lambda _self: None)

    sites_vcf = tmp_path / "sites.vcf"
    sites_vcf.write_text("##fileformat=VCFv4.2\n")
    fasta = tmp_path / "ref.fa"
    fasta.touch()
    chain = tmp_path / "chain.gz"
    chain.touch()

    # Set up a rejected.vcf with NoTarget so the dominant rejection is identified
    rej = tmp_path / "rejected.vcf"
    rej.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=249250621>\n"
        '##FILTER=<ID=NoTarget,Description="x">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "\n".join(
            f"chr1\t{i}\trs{i}\tA\tG\t.\tNoTarget\t."
            for i in range(1, 51)
        )
        + "\n"
    )

    with pytest.raises(LiftoverYieldError, match=r"50/100.*NoTarget"):
        lift_aadr_sites(
            sites_vcf_path=sites_vcf,
            chain_path=chain,
            target_fasta_path=fasta,
            output_lifted_vcf=tmp_path / "lifted.vcf",
            output_rejected_vcf=rej,
            input_filter_counters=Stage1InputFilters(
                palindrome_drops=0, non_snp_drops=0, non_autosome_drops=0, rows_written=100
            ),
        )


def test_lift_aadr_sites_threshold_warn_only(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """80% yield → above 70% fail threshold but below 95% warn threshold → warning logged."""
    import logging
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))

    import subprocess
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    fake_warn_stderr = """
    INFO	LiftoverVcf	Processed 100 variants.
    INFO	LiftoverVcf	20 variants failed to liftover.
    INFO	LiftoverVcf	0 variants lifted over but had mismatching reference alleles after lift over.
    INFO	LiftoverVcf	20.0000% of variants were not successfully lifted over
    INFO	LiftoverVcf	0 variants were lifted by swapping REF/ALT alleles.
    """

    from pileup_aadr import lift, tool_wrapper

    def fake_run(self: object, args: list[str], *, capture_stderr_to: Path, **_kw: object) -> object:
        capture_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        capture_stderr_to.write_text(fake_warn_stderr)
        return tool_wrapper.ToolRunResult(
            exit_code=0, stdout=None, stderr_path=capture_stderr_to,
            stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None,
        )

    monkeypatch.setattr(lift.ToolWrapper, "run", fake_run)
    monkeypatch.setattr(lift.ToolWrapper, "_check_version", lambda _self: None)
    caplog.set_level(logging.WARNING, logger="pileup_aadr.lift")

    sites_vcf = tmp_path / "sites.vcf"
    sites_vcf.write_text("##fileformat=VCFv4.2\n")
    fasta = tmp_path / "ref.fa"
    fasta.touch()
    chain = tmp_path / "chain.gz"
    chain.touch()
    rej = tmp_path / "rejected.vcf"  # leave nonexistent → empty bucket counts

    counters = lift_aadr_sites(
        sites_vcf_path=sites_vcf,
        chain_path=chain,
        target_fasta_path=fasta,
        output_lifted_vcf=tmp_path / "lifted.vcf",
        output_rejected_vcf=rej,
        input_filter_counters=Stage1InputFilters(
            palindrome_drops=0, non_snp_drops=0, non_autosome_drops=0, rows_written=100
        ),
    )
    assert counters.liftover_yield_pct == 80.0
    assert counters.liftover_yield_warning is True
    assert any("yield 80.00%" in r.message for r in caplog.records)
