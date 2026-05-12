"""Tests for pileup_call.py — Stage 3 (mpileup | pileupCaller pipe).

The Picard-equivalent integration tests mock `ToolWrapper.pipe` to return canned
exit codes + write a captured pileupCaller stderr fixture. Tests cover: clean
run → counters populated; downstream non-zero → ToolSubprocessError; upstream
non-zero (other than SIGPIPE-141) → ToolSubprocessError; upstream 141 +
downstream 0 → tolerated; stderr-parser format paths.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from pileup_aadr.errors import PileupAadrInternalError, ToolSubprocessError
from pileup_aadr.pileup_call import (
    parse_pileupcaller_stderr,
    run_pileup_call,
)
from pileup_aadr.tool_wrapper import ToolRunResult

FIXTURES_STDERR = Path(__file__).parent / "fixtures" / "stderr"


# --- parse_pileupcaller_stderr ---


def test_parse_clean_run() -> None:
    """Captured 1131600-site pileupCaller stderr parses to PileupCallerSummary."""
    text = (FIXTURES_STDERR / "pileupcaller_clean.stderr").read_text()
    summary = parse_pileupcaller_stderr(text)
    assert summary.total_sites == 1131600
    assert summary.non_missing_calls == 1107213
    assert summary.avg_raw_reads == 21.4
    assert summary.avg_damage_cleaned_reads == 21.4
    assert summary.avg_sampled_from == 21.4


def test_parse_missing_header_raises() -> None:
    """Stderr lacking the SampleName header → PileupAadrInternalError."""
    text = "Some unrelated output without the header.\n"
    with pytest.raises(PileupAadrInternalError, match="header"):
        parse_pileupcaller_stderr(text)


def test_parse_missing_data_line_raises() -> None:
    """Header present but no data line → PileupAadrInternalError."""
    text = (
        "SampleName\tTotalSites\tNonMissingCalls\tavgRawReads"
        "\tavgDamageCleanedReads\tavgSampledFrom\n"
        "# (no data row at all)\n"
    )
    with pytest.raises(PileupAadrInternalError, match="data line"):
        parse_pileupcaller_stderr(text)


# --- run_pileup_call (pipe mocked) ---


def _setup_run_mocks(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    *,
    upstream_exit: int,
    downstream_exit: int,
    pileupcaller_stderr_text: str,
    mpileup_stderr_text: str = "",
) -> None:
    """Patch subprocess.run (version probe) + ToolWrapper.pipe + _check_version."""
    import subprocess
    from unittest.mock import MagicMock

    # samtools/pileupCaller version probes both stub out
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="1.6.0.0", stderr="", returncode=0),
    )

    from pileup_aadr import pileup_call

    # Bypass binary-on-PATH lookup — neither samtools nor pileupCaller may be in
    # the test env. _resolve_binary returns a Path; the actual binary is never
    # invoked because we replace `pipe()` below.
    monkeypatch.setattr(
        pileup_call.ToolWrapper, "_resolve_binary",
        lambda _self, spec: Path(f"/usr/bin/fake_{spec.binary}"),
    )

    def fake_pipe(
        self: object,
        downstream: object,
        *,
        upstream_args: list[str],
        downstream_args: list[str],
        upstream_stderr_to: Path,
        downstream_stderr_to: Path,
    ) -> tuple[ToolRunResult, ToolRunResult]:
        upstream_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        upstream_stderr_to.write_text(mpileup_stderr_text)
        downstream_stderr_to.write_text(pileupcaller_stderr_text)
        return (
            ToolRunResult(
                exit_code=upstream_exit,
                stdout=None,
                stderr_path=upstream_stderr_to,
                stderr_text=None,
                wallclock_seconds=0.1,
                peak_rss_mb=None,
            ),
            ToolRunResult(
                exit_code=downstream_exit,
                stdout=None,
                stderr_path=downstream_stderr_to,
                stderr_text=None,
                wallclock_seconds=0.1,
                peak_rss_mb=None,
            ),
        )

    monkeypatch.setattr(pileup_call.ToolWrapper, "pipe", fake_pipe)
    monkeypatch.setattr(pileup_call.ToolWrapper, "_check_version", lambda _self: None)


def _common_run_args(tmp_path: Path) -> dict[str, Any]:
    bam = tmp_path / "user.bam"
    bam.touch()
    snp = tmp_path / "aadr.snp"
    snp.touch()
    bed = tmp_path / "aadr.bed"
    bed.touch()
    fasta = tmp_path / "ref.fa"
    fasta.touch()
    return {
        "bam_path": bam,
        "snp_path": snp,
        "bed_path": bed,
        "target_fasta_path": fasta,
        "output_prefix": tmp_path / "call" / "user_hg38",
        "sample_name": "Carsten",
        "pop_name": "TestPop",
    }


def test_clean_run_returns_counters(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Both processes exit 0 + pileupCaller stderr parses → Stage3CallCounters."""
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=0, downstream_exit=0,
        pileupcaller_stderr_text=(
            FIXTURES_STDERR / "pileupcaller_clean.stderr"
        ).read_text(),
    )
    counters = run_pileup_call(**_common_run_args(tmp_path))
    assert counters.pileupcaller_summary.total_sites == 1131600
    assert counters.pileupcaller_summary.non_missing_calls == 1107213
    assert counters.wallclock_seconds >= 0


def test_downstream_nonzero_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """pileupCaller exit != 0 → ToolSubprocessError naming pileupCaller."""
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=0, downstream_exit=1,
        pileupcaller_stderr_text="ERROR: bad input file\n",
    )
    with pytest.raises(ToolSubprocessError) as excinfo:
        run_pileup_call(**_common_run_args(tmp_path))
    assert "pileupCaller" in excinfo.value.what
    assert "exit code 1" in excinfo.value.why


def test_upstream_sigpipe_with_clean_downstream_tolerated(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Upstream exits 141 (SIGPIPE) + downstream exits 0 → tolerated, no raise."""
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=141, downstream_exit=0,
        pileupcaller_stderr_text=(
            FIXTURES_STDERR / "pileupcaller_clean.stderr"
        ).read_text(),
    )
    counters = run_pileup_call(**_common_run_args(tmp_path))
    assert counters.pileupcaller_summary.total_sites == 1131600


def test_upstream_real_error_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Upstream exit != 0 and != 141 → ToolSubprocessError naming samtools."""
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=1, downstream_exit=0,
        pileupcaller_stderr_text=(
            FIXTURES_STDERR / "pileupcaller_clean.stderr"
        ).read_text(),
        mpileup_stderr_text="samtools: missing BAM index\n",
    )
    with pytest.raises(ToolSubprocessError) as excinfo:
        run_pileup_call(**_common_run_args(tmp_path))
    assert "samtools mpileup" in excinfo.value.what
    assert "missing BAM index" in excinfo.value.why


def test_thread_cap_applied_by_default(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """threads=8 with cap → effective_threads becomes 4 (logged at INFO)."""
    import logging
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=0, downstream_exit=0,
        pileupcaller_stderr_text=(
            FIXTURES_STDERR / "pileupcaller_clean.stderr"
        ).read_text(),
    )
    caplog.set_level(logging.INFO, logger="pileup_aadr.pileup_call")
    run_pileup_call(threads=8, **_common_run_args(tmp_path))
    assert any("Capping mpileup threads 8 -> 4" in r.message for r in caplog.records)


def test_no_thread_cap_skips_capping(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """threads=8 with no_thread_cap=True → no cap log emitted."""
    import logging
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=0, downstream_exit=0,
        pileupcaller_stderr_text=(
            FIXTURES_STDERR / "pileupcaller_clean.stderr"
        ).read_text(),
    )
    caplog.set_level(logging.INFO, logger="pileup_aadr.pileup_call")
    run_pileup_call(threads=8, no_thread_cap=True, **_common_run_args(tmp_path))
    assert not any("Capping mpileup threads" in r.message for r in caplog.records)
