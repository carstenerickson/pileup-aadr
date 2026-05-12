"""Tests for concurrency.py — output_lock + tempdir + filesystem detection."""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from pileup_aadr.concurrency import (
    detect_filesystem_type,
    output_lock,
    tempdir,
    warn_if_networked_fs,
)
from pileup_aadr.errors import OutputLockHeldError

# --- output_lock ---


def test_output_lock_acquires_and_releases(tmp_path: Path) -> None:
    """Single-process acquire/release leaves a 0-byte .lock file behind."""
    prefix = tmp_path / "out"
    with output_lock(prefix):
        # Holder sidecar exists during the held window
        holder = Path(f"{prefix}.lock.holder")
        assert holder.exists()
        assert holder.read_text().strip() == str(os.getpid())
    # After release: lock file persists (0-byte), holder removed
    assert Path(f"{prefix}.lock").exists()
    assert not Path(f"{prefix}.lock.holder").exists()


def test_output_lock_creates_parent_dir(tmp_path: Path) -> None:
    """Lock acquisition mkdir-p's the parent dir."""
    prefix = tmp_path / "deep" / "nested" / "out"
    with output_lock(prefix):
        pass
    assert Path(f"{prefix}.lock").exists()


def test_output_lock_contention_raises(tmp_path: Path) -> None:
    """Second acquire while first held → OutputLockHeldError with PID diagnostic.

    Uses a real subprocess.Popen to hold the lock. The subprocess writes a
    'ready' sentinel file once it's holding the lock, then sleeps. We poll
    for the sentinel, then attempt to acquire and verify the error names the
    subprocess PID.
    """
    prefix = tmp_path / "out"
    sentinel = tmp_path / "ready.flag"
    holder_script = (
        "import sys, time, fcntl, os\n"
        "from pathlib import Path\n"
        f"sys.path.insert(0, {str(Path('.').resolve())!r})\n"
        "from pileup_aadr.concurrency import output_lock\n"
        f"with output_lock(Path({str(prefix)!r})):\n"
        f"    Path({str(sentinel)!r}).write_text('ready')\n"
        "    time.sleep(3.0)\n"
    )
    proc = subprocess.Popen([sys.executable, "-c", holder_script])
    try:
        # Wait for the holder to acquire
        deadline = time.monotonic() + 5
        while not sentinel.exists() and time.monotonic() < deadline:
            time.sleep(0.05)
        assert sentinel.exists(), "subprocess never signaled lock acquisition"

        with pytest.raises(OutputLockHeldError) as excinfo:
            with output_lock(prefix):
                pass
        assert "PID" in excinfo.value.why
        assert str(proc.pid) in excinfo.value.why
    finally:
        proc.terminate()
        proc.wait(timeout=3)


# --- tempdir ---


def test_tempdir_clean_exit_removes_dir(tmp_path: Path) -> None:
    """Default tempdir context: dir is removed on clean exit."""
    captured: list[Path] = []
    with tempdir(base=tmp_path) as td:
        captured.append(td)
        assert td.exists()
    assert not captured[0].exists()


def test_tempdir_keeps_on_crash_by_default(tmp_path: Path) -> None:
    """Default behavior: any uncaught exception RETAINS the tempdir for forensics."""
    captured: list[Path] = []
    with pytest.raises(RuntimeError, match="boom"):
        with tempdir(base=tmp_path) as td:
            captured.append(td)
            raise RuntimeError("boom")
    assert captured[0].exists()  # retained for forensics


def test_tempdir_keep_always_retains_on_clean_exit(tmp_path: Path) -> None:
    """keep_always=True: dir survives clean exit."""
    captured: list[Path] = []
    with tempdir(base=tmp_path, keep_always=True) as td:
        captured.append(td)
    assert captured[0].exists()


def test_tempdir_clean_on_crash_removes_dir(tmp_path: Path) -> None:
    """clean_on_crash=True: even crash path removes the dir."""
    captured: list[Path] = []
    with pytest.raises(RuntimeError, match="boom"):
        with tempdir(base=tmp_path, clean_on_crash=True) as td:
            captured.append(td)
            raise RuntimeError("boom")
    assert not captured[0].exists()


def test_tempdir_uses_pileup_aadr_prefix(tmp_path: Path) -> None:
    """Tempdir name starts with `pileup-aadr-` for identifiability in `ls /tmp`."""
    with tempdir(base=tmp_path) as td:
        assert td.name.startswith("pileup-aadr-")


# --- filesystem-type detection ---


def test_detect_filesystem_type_returns_none_on_non_linux(tmp_path: Path) -> None:
    """Non-Linux: detect_filesystem_type returns None (not implemented)."""
    if sys.platform == "linux":
        pytest.skip("Linux: detect_filesystem_type returns the actual fstype")
    assert detect_filesystem_type(tmp_path) is None


def test_warn_if_networked_fs_no_warn_on_local(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Local-FS path (or non-Linux): no NFS/SMB warning emitted."""
    import logging
    caplog.set_level(logging.WARNING, logger="pileup_aadr.concurrency")
    warn_if_networked_fs(tmp_path / "out")
    assert not any("nfs" in r.message.lower() for r in caplog.records)
