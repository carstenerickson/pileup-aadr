"""Tests for dict_resolve.py — target FASTA `.dict` lookup + auto-generation."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from pileup_aadr.dict_resolve import (
    ensure_target_fasta_dict,
    find_existing_dict,
    find_or_user_cache_dict_path,
)
from pileup_aadr.tool_wrapper import ToolRunResult

# --- find_existing_dict ---


def test_find_existing_long_form(tmp_path: Path) -> None:
    """`<fasta>.dict` (Picard's default output) is preferred when present."""
    fasta = tmp_path / "hg38.fa"
    fasta.touch()
    long_form = Path(f"{fasta}.dict")
    long_form.touch()
    assert find_existing_dict(fasta) == long_form


def test_find_existing_short_form(tmp_path: Path) -> None:
    """`<fasta-stem>.dict` (GATK resource bundle convention) is the fallback."""
    fasta = tmp_path / "hg38.fa"
    fasta.touch()
    short_form = fasta.with_suffix(".dict")
    short_form.touch()
    assert find_existing_dict(fasta) == short_form


def test_find_existing_returns_none_when_neither_present(tmp_path: Path) -> None:
    fasta = tmp_path / "hg38.fa"
    fasta.touch()
    assert find_existing_dict(fasta) is None


def test_find_existing_long_form_wins_when_both_present(tmp_path: Path) -> None:
    """If both candidates exist, prefer the long form."""
    fasta = tmp_path / "hg38.fa"
    fasta.touch()
    long_form = Path(f"{fasta}.dict")
    short_form = fasta.with_suffix(".dict")
    long_form.touch()
    short_form.touch()
    assert find_existing_dict(fasta) == long_form


# --- find_or_user_cache_dict_path ---


def test_cache_path_alongside_when_writable(tmp_path: Path) -> None:
    """Writable parent dir → `.dict` goes alongside the FASTA."""
    fasta = tmp_path / "hg38.fa"
    fasta.touch()
    out = find_or_user_cache_dict_path(fasta)
    assert out == Path(f"{fasta.resolve()}.dict")


def test_cache_path_user_cache_when_readonly(tmp_path: Path) -> None:
    """Read-only parent → cache under the user-cache dir, keyed by SHA of abs path."""
    fasta = tmp_path / "ro_dir" / "hg38.fa"
    fasta.parent.mkdir()
    fasta.touch()
    fasta.parent.chmod(0o555)  # read+execute only
    cache = tmp_path / "test_cache"
    try:
        out = find_or_user_cache_dict_path(fasta, cache_dir=cache)
        assert out.parent == cache
        assert out.name.startswith("hg38.fa.")
        assert out.name.endswith(".dict")
        assert cache.exists()
    finally:
        fasta.parent.chmod(0o755)  # restore for cleanup


# --- ensure_target_fasta_dict ---


def test_ensure_returns_existing_without_invoking_picard(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """When `.dict` exists, no Picard subprocess is invoked."""
    fasta = tmp_path / "hg38.fa"
    fasta.touch()
    long_form = Path(f"{fasta}.dict")
    long_form.touch()

    from pileup_aadr import dict_resolve

    # Sentinel — if ToolWrapper is constructed at all, the test fails
    def boom(*_a: object, **_kw: object) -> object:
        raise AssertionError("ToolWrapper should not be invoked when .dict exists")
    monkeypatch.setattr(dict_resolve, "ToolWrapper", boom)

    assert ensure_target_fasta_dict(fasta) == long_form


def test_ensure_invokes_picard_when_dict_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """Missing `.dict` triggers Picard CreateSequenceDictionary; output path returned."""
    fasta = tmp_path / "hg38.fa"
    fasta.touch()

    import subprocess
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="Version:3.3.0", stderr="", returncode=0),
    )

    from pileup_aadr import dict_resolve

    monkeypatch.setattr(
        dict_resolve.ToolWrapper, "_resolve_binary",
        lambda _self, _spec: tmp_path / "fake_picard.jar",
    )
    monkeypatch.setattr(
        dict_resolve.ToolWrapper, "_check_version", lambda _self: None,
    )

    captured: dict[str, list[str]] = {}

    def fake_run(
        self: object,
        *,
        args: list[str],
        capture_stderr_to: Path,
        check: bool = False,
        **_kw: object,
    ) -> ToolRunResult:
        captured["args"] = args
        # Simulate Picard writing the .dict
        out_idx = args.index("--OUTPUT") + 1
        Path(args[out_idx]).write_text("@HD\tVN:1.6\n")
        capture_stderr_to.write_text("")
        return ToolRunResult(
            exit_code=0, stdout=None, stderr_path=capture_stderr_to,
            stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None,
        )

    monkeypatch.setattr(dict_resolve.ToolWrapper, "run", fake_run)

    result = ensure_target_fasta_dict(fasta)
    assert result == Path(f"{fasta.resolve()}.dict")
    assert result.exists()
    # Confirm Picard was invoked with the expected subcommand + args
    assert captured["args"][:3] == ["CreateSequenceDictionary", "--REFERENCE", str(fasta)]
