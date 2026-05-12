"""Tests for tool_wrapper.py — ToolSpec / ToolRunResult / ToolWrapper.

Heavy reliance on mocking subprocess.run since real binary invocation is gated
behind @pytest.mark.requires_* markers (Day 2 doesn't have a fixture BAM yet).
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from pileup_aadr.errors import (
    JavaNotFoundError,
    PicardNotFoundError,
    PileupCallerNotFoundError,
    SamtoolsNotFoundError,
    ToolVersionError,
)
from pileup_aadr.tool_wrapper import (
    JAVA_SPEC,
    MOSDEPTH_SPEC,
    PICARD_SPEC,
    PILEUPCALLER_SPEC,
    SAMTOOLS_SPEC,
    ToolRunResult,
    ToolSpec,
    ToolWrapper,
)


def _clear_version_cache() -> None:
    """Clear the class-level version cache between tests (mocks don't share)."""
    ToolWrapper._version_cache.clear()


@pytest.fixture(autouse=True)
def _reset_cache() -> None:
    _clear_version_cache()


# --- ToolSpec ---


def test_all_5_specs_have_required_fields() -> None:
    for spec in (SAMTOOLS_SPEC, PILEUPCALLER_SPEC, PICARD_SPEC, JAVA_SPEC, MOSDEPTH_SPEC):
        assert spec.binary
        assert spec.version_args
        assert spec.version_regex
        assert spec.min_version
        assert spec.tested_against


def test_picard_spec_is_jar() -> None:
    assert PICARD_SPEC.is_jar is True


def test_non_jar_specs_default_is_jar_false() -> None:
    for spec in (SAMTOOLS_SPEC, PILEUPCALLER_SPEC, JAVA_SPEC, MOSDEPTH_SPEC):
        assert spec.is_jar is False


# --- ToolWrapper binary resolution + version probe (mocked subprocess) ---


def test_samtools_not_on_path_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda _: None)
    with pytest.raises(SamtoolsNotFoundError, match="not found on PATH"):
        ToolWrapper(SAMTOOLS_SPEC)


def test_pileupcaller_not_on_path_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda _: None)
    with pytest.raises(PileupCallerNotFoundError):
        ToolWrapper(PILEUPCALLER_SPEC)


def test_java_not_on_path_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda _: None)
    with pytest.raises(JavaNotFoundError):
        ToolWrapper(JAVA_SPEC)


def test_picard_jar_not_findable_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """No $PICARD_JAR, no $CONDA_PREFIX, no wrapper, no known paths → PicardNotFoundError."""
    monkeypatch.delenv("PICARD_JAR", raising=False)
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    monkeypatch.setattr("shutil.which", lambda _: None)
    monkeypatch.setattr(Path, "exists", lambda _: False)
    with pytest.raises(PicardNotFoundError, match="not found"):
        ToolWrapper(PICARD_SPEC, skip_version_check=True)


def test_picard_jar_via_env_var(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """$PICARD_JAR pointing at a real .jar resolves cleanly."""
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"PK\x03\x04 fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))
    # Mock subprocess.run to return a fake version string
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )
    wrapper = ToolWrapper(PICARD_SPEC)
    assert wrapper.binary_path == fake_jar
    assert wrapper.version() == "3.3.0"


def test_version_too_old_raises(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Picard 2.27.5 < 3.0.0 → ToolVersionError."""
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:2.27.5\n", returncode=0),
    )
    with pytest.raises(ToolVersionError, match="version too old"):
        ToolWrapper(PICARD_SPEC)


def test_version_unparseable_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Tool emits no parseable version string → ToolVersionError with format-changed hint."""
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="garbage", stderr="more garbage", returncode=0),
    )
    with pytest.raises(ToolVersionError, match="matched no line"):
        ToolWrapper(PICARD_SPEC)


def test_version_cached_per_binary(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Second call to .version() returns cached value (no second subprocess invocation)."""
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))

    call_count = 0

    def counting_run(*_args: object, **_kwargs: object) -> object:
        nonlocal call_count
        call_count += 1
        return MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0)

    monkeypatch.setattr(subprocess, "run", counting_run)

    wrapper = ToolWrapper(PICARD_SPEC, skip_version_check=True)
    assert wrapper.version() == "3.3.0"
    assert wrapper.version() == "3.3.0"  # cached
    # Cached: only ONE subprocess invocation across two version() calls
    assert call_count == 1


# --- _build_invocation ---


def test_build_invocation_jar_prepends_java_jar(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """is_jar specs build `["java", *jvm_args, "-jar", <jar>, *args]`."""
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    wrapper = ToolWrapper(PICARD_SPEC, skip_version_check=True)
    cmd = wrapper._build_invocation(["LiftoverVcf", "--CHAIN", "/tmp/c.gz"])
    assert cmd == ["java", "-jar", str(fake_jar), "LiftoverVcf", "--CHAIN", "/tmp/c.gz"]


def test_build_invocation_jar_with_jvm_args(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """jvm_args injected between `java` and `-jar` (Day-2 H9 fix)."""
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    wrapper = ToolWrapper(PICARD_SPEC, skip_version_check=True)
    cmd = wrapper._build_invocation(
        ["LiftoverVcf", "--CHAIN", "/tmp/c.gz"],
        jvm_args=["-Xmx3g"],
    )
    assert cmd == ["java", "-Xmx3g", "-jar", str(fake_jar), "LiftoverVcf", "--CHAIN", "/tmp/c.gz"]


def test_build_invocation_non_jar_with_jvm_args_raises(
    monkeypatch: pytest.MonkeyPatch
) -> None:
    """Non-jar specs reject jvm_args (defensive)."""
    monkeypatch.setattr("shutil.which", lambda _b: f"/usr/local/bin/{_b}")
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="samtools 1.23.1\n", stderr="", returncode=0),
    )

    wrapper = ToolWrapper(SAMTOOLS_SPEC, skip_version_check=True)
    with pytest.raises(ValueError, match="jvm_args is only valid for is_jar"):
        wrapper._build_invocation(["mpileup"], jvm_args=["-Xmx3g"])


# --- ToolRunResult dataclass ---


def test_toolrunresult_is_frozen() -> None:
    """ToolRunResult is frozen — caller mutation should raise FrozenInstanceError."""
    result = ToolRunResult(
        exit_code=0,
        stdout=None,
        stderr_path=Path("/tmp/x"),
        stderr_text=None,
        wallclock_seconds=1.0,
        peak_rss_mb=None,
    )
    with pytest.raises(Exception):  # noqa: B017 — FrozenInstanceError is dataclasses-private
        result.exit_code = 1  # type: ignore[misc]


# --- ToolSpec dataclass ---


def test_toolspec_is_frozen() -> None:
    spec = ToolSpec(
        binary="x",
        version_args=["-v"],
        version_regex=r"(\d+)",
        min_version="1",
        tested_against="1",
        error_class_missing=SamtoolsNotFoundError,
    )
    with pytest.raises(Exception):  # noqa: B017
        spec.binary = "y"  # type: ignore[misc]
