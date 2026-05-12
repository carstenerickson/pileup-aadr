"""Tests for `extract_cmd` — the click wrapper around `run_extract`.

Mostly thin contract tests for the CLI <-> orchestrator boundary. Catches the
class of bug where a CLI flag's semantics drift from the orchestrator's
internal field representation (the BAQ-flip inversion was the motivating
example — see test_enable_baq_flag_maps_correctly_to_orchestrator).
"""
from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from pileup_aadr.cli import cli
from pileup_aadr.types import ExtractCliArgs


def _make_inputs(tmp_path: Path) -> tuple[Path, Path]:
    """Touch the two positional inputs so click's `exists=True` path checks pass."""
    bam = tmp_path / "user.bam"
    bam.touch()
    aadr = tmp_path / "aadr.snp"
    aadr.touch()
    return bam, aadr


def _captured_args(monkeypatch: pytest.MonkeyPatch) -> dict[str, ExtractCliArgs]:
    """Patch run_extract to capture the ExtractCliArgs without executing."""
    captured: dict[str, ExtractCliArgs] = {}

    def fake_run(args: ExtractCliArgs) -> int:
        captured["args"] = args
        return 0

    from pileup_aadr import extract_cmd
    monkeypatch.setattr(extract_cmd, "run_extract", fake_run)
    return captured


# --- BAQ flag mapping (regression for the v0.2 LLD #19 inversion bug) ---


def test_enable_baq_flag_default_disables_baq_in_mpileup(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """CLI default (no --enable-baq) must produce no_baq=False so the
    orchestrator appends `-B` to mpileup args, DISABLING samtools BAQ —
    which is the HLD-spec'd default ("default: -B is passed, disabling
    samtools BAQ to match pileupCaller's recommended cmdline").

    Inversion regression: the v0.1.0-0.1.2 CLI used `not enable_baq` which
    produced no_baq=True by default → -B NOT appended → BAQ ENABLED
    (opposite of HLD). The bug was caught by LLD #19 layer-B and only
    surfaced when the CLI path was compared against direct ExtractCliArgs
    instantiation (which had no_baq=False as the dataclass default).
    """
    bam, aadr = _make_inputs(tmp_path)
    captured = _captured_args(monkeypatch)
    runner = CliRunner()

    result = runner.invoke(
        cli, ["extract", str(bam), str(aadr), "-o", str(tmp_path / "out")],
    )
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    assert captured["args"].no_baq is False, (
        f"CLI default should produce no_baq=False (HLD: -B passed by default, "
        f"BAQ disabled) — got no_baq={captured['args'].no_baq}"
    )


def test_enable_baq_flag_set_enables_baq_in_mpileup(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """`--enable-baq` must produce no_baq=True so the orchestrator does
    NOT append `-B`, leaving samtools BAQ ENABLED (the opt-in case)."""
    bam, aadr = _make_inputs(tmp_path)
    captured = _captured_args(monkeypatch)
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "extract", str(bam), str(aadr),
            "-o", str(tmp_path / "out"), "--enable-baq",
        ],
    )
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    assert captured["args"].no_baq is True, (
        f"--enable-baq must produce no_baq=True (BAQ enabled in mpileup) — "
        f"got no_baq={captured['args'].no_baq}"
    )


def test_dataclass_default_matches_cli_default(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """The dataclass default and the CLI no-flag default must produce the
    same `no_baq` value. Catches the class of bug where the CLI layer's
    flip + the dataclass default get inverted relative to each other (which
    the LLD #19 layer-B test surfaced as a Stage-3 .geno divergence).
    """
    bam, aadr = _make_inputs(tmp_path)
    captured = _captured_args(monkeypatch)
    runner = CliRunner()

    runner.invoke(
        cli, ["extract", str(bam), str(aadr), "-o", str(tmp_path / "out")],
    )
    cli_default_no_baq = captured["args"].no_baq

    # Dataclass default — what direct programmatic instantiation produces
    dataclass_default = ExtractCliArgs(
        bam=bam, aadr_snp=aadr, output_prefix=tmp_path / "out",
    )
    assert cli_default_no_baq == dataclass_default.no_baq, (
        f"CLI no-flag default ({cli_default_no_baq}) must equal dataclass "
        f"default ({dataclass_default.no_baq}) — divergence here means "
        f"`pileup-aadr extract` and direct `run_extract(ExtractCliArgs(...))` "
        f"produce different mpileup args, which produces different .geno "
        f"output, which produces different f2 numbers downstream."
    )
