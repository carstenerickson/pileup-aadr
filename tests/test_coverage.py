"""Tests for `coverage` subcommand — mosdepth wrapper + summary parser."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import pytest
from click.testing import CliRunner

from pileup_aadr.cli import cli
from pileup_aadr.coverage_impl import _parse_mosdepth_summary, run_coverage
from pileup_aadr.tool_wrapper import ToolRunResult
from pileup_aadr.types import CoverageCliArgs

# --- _parse_mosdepth_summary ---


def test_parse_summary_basic_per_chrom(tmp_path: Path) -> None:
    """Round-trip a 3-row mosdepth summary into the expected dict shape."""
    summary = tmp_path / "out.mosdepth.summary.txt"
    summary.write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        "chr1\t249250621\t12345678\t30.5\t0\t100\n"
        "chr22\t51304566\t2000000\t25.0\t0\t80\n"
        "total\t300555187\t14345678\t29.6\t0\t100\n"
    )
    parsed = _parse_mosdepth_summary(summary)
    assert parsed["per_chrom"]["chr1"] == {
        "length": 249_250_621, "bases": 12_345_678, "mean_coverage": 30.5,
    }
    assert parsed["per_chrom"]["chr22"]["mean_coverage"] == 25.0
    assert parsed["per_chrom"]["total"]["bases"] == 14_345_678


def test_parse_summary_skips_short_rows(tmp_path: Path) -> None:
    """Defensive: rows with < 6 cols (corruption) are skipped, not crash."""
    summary = tmp_path / "out.mosdepth.summary.txt"
    summary.write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        "chr1\t249250621\t12345678\t30.5\t0\t100\n"
        "truncated_row_only_3_cols\t10\t20\n"
    )
    parsed = _parse_mosdepth_summary(summary)
    assert "chr1" in parsed["per_chrom"]
    assert "truncated_row_only_3_cols" not in parsed["per_chrom"]


# --- run_coverage (mocked mosdepth subprocess) ---


def _setup_mosdepth_mock(
    monkeypatch: pytest.MonkeyPatch,
    *,
    summary_text: str,
    dist_text: str = "",
) -> None:
    """Patch ToolWrapper to bypass binary lookup + simulate mosdepth's outputs.

    Both summary.txt and global.dist.txt are written under <prefix>.mosdepth.*
    in the orchestrator's tempdir so coverage_impl's parsers can read them.
    """
    import subprocess
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="mosdepth 0.3.6", stderr="", returncode=0),
    )

    from pileup_aadr import coverage_impl

    monkeypatch.setattr(
        coverage_impl.ToolWrapper, "_resolve_binary",
        lambda _self, spec: Path(f"/usr/bin/fake_{spec.binary}"),
    )
    monkeypatch.setattr(
        coverage_impl.ToolWrapper, "_check_version", lambda _self: None,
    )

    def fake_run(
        self: object,
        *,
        args: list[str],
        capture_stderr_to: Path,
        check: bool = False,
        **_kw: Any,
    ) -> ToolRunResult:
        # Find the prefix mosdepth would write to (last positional before BAM)
        # The args list ends with [..., prefix, bam] per coverage_impl
        prefix = Path(args[-2])
        summary_path = Path(f"{prefix}.mosdepth.summary.txt")
        dist_path = Path(f"{prefix}.mosdepth.global.dist.txt")
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(summary_text)
        dist_path.write_text(dist_text)
        capture_stderr_to.write_text("")
        return ToolRunResult(
            exit_code=0, stdout=None, stderr_path=capture_stderr_to,
            stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None,
        )

    monkeypatch.setattr(coverage_impl.ToolWrapper, "run", fake_run)


# mosdepth global.dist.txt format: chrom\tdepth\tcumulative_fraction (descending depth).
# Synthetic distribution: 30x covers 50% of chr1, 10x covers 80%, 5x covers 95%, 1x = 99%.
_DIST_CHR1 = (
    "chr1\t100\t0.0\n"
    "chr1\t30\t0.5\n"
    "chr1\t10\t0.8\n"
    "chr1\t5\t0.95\n"
    "chr1\t1\t0.99\n"
    "chr1\t0\t1.0\n"
)
_DIST_CHR22 = (
    "chr22\t30\t0.1\n"
    "chr22\t10\t0.4\n"
    "chr22\t5\t0.7\n"
    "chr22\t1\t0.95\n"
    "chr22\t0\t1.0\n"
)


def test_run_coverage_emits_tsv(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Default (TSV) mode: 9-col header + per-chrom row with median + frac>=Nx."""
    bam = tmp_path / "user.bam"
    bam.touch()
    _setup_mosdepth_mock(
        monkeypatch,
        summary_text=(
            "chrom\tlength\tbases\tmean\tmin\tmax\n"
            "chr1\t249250621\t12345678\t30.5\t0\t100\n"
            "chr22\t51304566\t2000000\t25.0\t0\t80\n"
        ),
        dist_text=_DIST_CHR1 + _DIST_CHR22,
    )
    args = CoverageCliArgs(bam=bam)
    exit_code = run_coverage(args)
    assert exit_code == 0

    out = capsys.readouterr().out.splitlines()
    # Header lists all 9 columns per HLD §"CLI reference > coverage"
    assert out[0] == (
        "chrom\tlength\tbases\tmean_coverage\tmedian_coverage"
        "\tfraction_at_>=1x\tfraction_at_>=5x"
        "\tfraction_at_>=10x\tfraction_at_>=30x"
    )
    # chr1: median is 30 (largest depth where frac>=0.5);
    # frac@1x=0.99, @5x=0.95, @10x=0.8, @30x=0.5
    assert any(
        line.startswith("chr1\t249250621\t12345678\t30.5\t30\t0.99\t0.95\t0.8\t0.5")
        for line in out
    )
    # chr22: median is 5 (largest depth where frac>=0.5; 0.7 at 5x, 0.4 at 10x)
    assert any(
        line.startswith("chr22\t51304566\t2000000\t25.0\t5\t0.95\t0.7\t0.4\t0.1")
        for line in out
    )


def test_run_coverage_emits_json(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """--json mode: per-chrom dict carries the median + 4 fraction fields."""
    bam = tmp_path / "user.bam"
    bam.touch()
    _setup_mosdepth_mock(
        monkeypatch,
        summary_text=(
            "chrom\tlength\tbases\tmean\tmin\tmax\n"
            "chr1\t249250621\t12345678\t30.5\t0\t100\n"
        ),
        dist_text=_DIST_CHR1,
    )
    args = CoverageCliArgs(bam=bam, json_output=True)
    run_coverage(args)
    payload = json.loads(capsys.readouterr().out)
    chr1 = payload["per_chrom"]["chr1"]
    assert chr1["mean_coverage"] == 30.5
    assert chr1["median_coverage"] == 30
    assert chr1["fraction_at_>=1x"] == 0.99
    assert chr1["fraction_at_>=5x"] == 0.95
    assert chr1["fraction_at_>=10x"] == 0.8
    assert chr1["fraction_at_>=30x"] == 0.5


def test_empty_global_dist_yields_nan_quantile_cols(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Defensive: empty global.dist.txt → median + fraction cols stringify
    as 'nan', not crash."""
    bam = tmp_path / "user.bam"
    bam.touch()
    _setup_mosdepth_mock(
        monkeypatch,
        summary_text=(
            "chrom\tlength\tbases\tmean\tmin\tmax\n"
            "chr1\t249250621\t12345678\t30.5\t0\t100\n"
        ),
        dist_text="",
    )
    args = CoverageCliArgs(bam=bam)
    run_coverage(args)
    out = capsys.readouterr().out.splitlines()
    # When dist is empty, _merge sets NaN; format() renders NaN as "nan" and
    # repr-floats it (median 0 since the for-loop never matches frac>=0.5).
    # Median goes through `int = 0`; the four frac cols are NaN.
    chr1_row = next(line for line in out if line.startswith("chr1\t"))
    parts = chr1_row.split("\t")
    # cols: chrom length bases mean median frac>=1 frac>=5 frac>=10 frac>=30
    # All five derived cols are NaN when there's no global.dist data
    for v in parts[4:9]:
        assert v == "nan"


# --- cli.py wiring ---


def test_cli_coverage_help_lists_command() -> None:
    """`pileup-aadr --help` mentions the coverage subcommand."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "coverage" in result.output


def test_cli_coverage_help_renders_options() -> None:
    """`pileup-aadr coverage --help` shows --regions / --threads / --json."""
    runner = CliRunner()
    result = runner.invoke(cli, ["coverage", "--help"])
    assert result.exit_code == 0
    assert "--regions" in result.output
    assert "--threads" in result.output
    assert "--json" in result.output
