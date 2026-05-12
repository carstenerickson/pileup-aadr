"""Concurrency primitives: output-prefix lock + tempdir lifecycle.

Two context managers used by `extract_orch.run_extract`:

- `output_lock` — advisory exclusive flock on `<prefix>.lock`; holder PID
  written to `<prefix>.lock.holder` for diagnostics. v2.1 critique H8 fix:
  PID lives in a separate sidecar file so we don't truncate the lock file
  itself on release (which masked the diagnostic for the next contender).

- `tempdir` — per-invocation tempdir with crash-survival semantics. Default:
  retain on crash for forensics, clean on success. `keep_always` and
  `clean_on_crash` flags override both directions.

`warn_if_networked_fs` emits a stderr WARNING when the output prefix is on
NFS/SMB where flock semantics may be silently no-op'd.
"""
from __future__ import annotations

import fcntl
import logging
import os
import shutil
import sys
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from .errors import OutputLockHeldError

log = logging.getLogger(__name__)


@contextmanager
def output_lock(prefix: Path) -> Iterator[Path]:
    """Acquire an advisory exclusive lock on `<prefix>.lock`.

    The lock file is created if missing; left in place on context exit
    (release is via fd close, not file unlink — leaving a 0-byte `.lock`
    file is harmless and avoids the classic unlink race).

    PID diagnostic: a separate `<prefix>.lock.holder` sidecar is written +
    removed alongside lock acquisition. Reading it can race (holder may
    have just released), but the diagnostic is best-effort — the real
    semantics rely on `flock`, not file content.

    Raises:
        OutputLockHeldError: another process holds the lock.
    """
    lock_path = Path(f"{prefix}.lock")
    holder_path = Path(f"{prefix}.lock.holder")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as e:
            holder_pid = _read_holder_pid(holder_path)
            holder_str = f"PID {holder_pid}" if holder_pid else "another process"
            raise OutputLockHeldError(
                what=str(lock_path),
                why=f"advisory lock held by {holder_str}",
                fix=(
                    f"Wait for {holder_str} to finish, or choose a different -o prefix"
                ),
            ) from e
        try:
            holder_path.write_text(f"{os.getpid()}\n")
        except OSError:
            log.debug("Could not write holder sidecar at %s (non-fatal)", holder_path)
        log.debug("Acquired output lock: %s", lock_path)
        yield lock_path
    finally:
        try:
            holder_path.unlink(missing_ok=True)
        except OSError:
            pass
        os.close(fd)
        log.debug("Released output lock: %s", lock_path)


def _read_holder_pid(holder_path: Path) -> str | None:
    """Best-effort read of the lock holder's PID. Returns None on any failure."""
    try:
        content = holder_path.read_text().strip()
        return content if content else None
    except OSError:
        return None


@contextmanager
def tempdir(
    *,
    base: Path | None = None,
    keep_always: bool = False,
    clean_on_crash: bool = False,
) -> Iterator[Path]:
    """Per-invocation tempdir with crash-survival semantics.

    Default behavior (neither flag set):
        - Clean exit: tempdir removed
        - Crash (any uncaught exception): tempdir RETAINED for forensics;
          path logged at ERROR.
    """
    base_dir = base or Path(tempfile.gettempdir())
    base_dir.mkdir(parents=True, exist_ok=True)
    td = Path(tempfile.mkdtemp(prefix="pileup-aadr-", dir=base_dir))
    log.info("Tempdir: %s", td)
    try:
        yield td
    except BaseException:
        if keep_always or not clean_on_crash:
            log.error("Tempdir RETAINED on crash for forensics: %s", td)
        else:
            shutil.rmtree(td, ignore_errors=True)
            log.warning("Tempdir cleaned on crash (--clean-tempdir-on-crash): %s", td)
        raise
    else:
        if keep_always:
            log.info("Tempdir RETAINED (--keep-tempdir): %s", td)
        else:
            shutil.rmtree(td, ignore_errors=True)
            log.debug("Tempdir cleaned: %s", td)


def detect_filesystem_type(path: Path) -> str | None:
    """Best-effort filesystem-type detection for the output prefix's parent dir.

    Used at startup to emit a stderr WARNING if output is on NFS/SMB (where
    `fcntl.flock` semantics may be silently no-op'd, breaking `output_lock`'s
    guarantee). Linux-only; macOS returns None.
    """
    if sys.platform != "linux":
        return None
    target = path.resolve()
    try:
        with open("/proc/mounts") as f:
            mounts = [line.split()[:3] for line in f if line.strip()]
    except OSError:
        return None
    mounts.sort(key=lambda m: len(m[1]), reverse=True)
    for _device, mount_point, fstype in mounts:
        if str(target).startswith(mount_point):
            return fstype
    return None


def warn_if_networked_fs(output_prefix: Path) -> None:
    """Emit stderr WARNING if output is on NFS/SMB/CIFS where flock may be no-op."""
    fstype = detect_filesystem_type(output_prefix.parent)
    if fstype and fstype in ("nfs", "nfs4", "cifs", "smbfs"):
        log.warning(
            "Output prefix is on %s filesystem; fcntl.flock semantics may be "
            "unreliable. For cluster workloads, wrap with workload-manager "
            "coordination (Slurm jobid, lockdir on local-fs scratch).",
            fstype,
        )


__all__ = [
    "detect_filesystem_type",
    "output_lock",
    "tempdir",
    "warn_if_networked_fs",
]
