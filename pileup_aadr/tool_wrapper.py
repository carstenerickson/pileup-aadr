"""Subprocess plumbing for the external tools (Picard, samtools, pileupCaller, mosdepth, java).

All three Stage modules + the `coverage` subcommand flow through ToolWrapper to get
uniform error handling, version checking, RSS measurement, and stderr-to-disk capture
per the HLD §"Subprocess wrapper abstraction" + §"Logging architecture" disciplines.
"""
from __future__ import annotations

import glob
import logging
import os
import re
import resource
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, ClassVar

from packaging.version import InvalidVersion, Version

from .errors import (
    JavaNotFoundError,
    MosdepthNotFoundError,
    PicardNotFoundError,
    PileupCallerNotFoundError,
    SamtoolsNotFoundError,
    ToolSubprocessError,
    ToolVersionError,
)

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class ToolSpec:
    """Per-tool metadata for ToolWrapper."""

    binary: str
    version_args: list[str]
    version_regex: str
    min_version: str
    tested_against: str
    error_class_missing: type[Exception]
    error_class_version: type[Exception] = field(default=ToolVersionError)
    is_jar: bool = False


@dataclass(frozen=True)
class ToolRunResult:
    """Result from ToolWrapper.run."""

    exit_code: int
    stdout: str | None  # only populated if explicitly captured
    stderr_path: Path  # always written; tail-readable for diagnostic
    stderr_text: str | None  # populated only when explicitly read
    wallclock_seconds: float
    peak_rss_mb: float | None  # None if measurement unavailable


# --- Tool specs (versions verified empirically v2.1) ---

SAMTOOLS_SPEC = ToolSpec(
    binary="samtools",
    version_args=["--version"],
    version_regex=r"^samtools\s+(\d+\.\d+(?:\.\d+)?)",
    min_version="1.16",
    tested_against="1.23.1",
    error_class_missing=SamtoolsNotFoundError,
)

PILEUPCALLER_SPEC = ToolSpec(
    binary="pileupCaller",
    version_args=["--version"],
    version_regex=r"(\d+\.\d+\.\d+\.\d+)",  # matches "1.6.0.0"
    min_version="1.6.0.0",
    tested_against="1.6.0.0",
    error_class_missing=PileupCallerNotFoundError,
)

PICARD_SPEC = ToolSpec(
    binary="picard.jar",  # JAR path resolved at startup; see resolve_picard_jar()
    # `java -jar picard.jar --version` (no subcommand) prints help + exits 1.
    # Version is only printed via a subcommand: `LiftoverVcf --version` -> "Version:3.1.1"
    # on stderr. We probe via LiftoverVcf since that's the one we'll actually invoke.
    version_args=["LiftoverVcf", "--version"],
    version_regex=r"Version:(\d+\.\d+(?:\.\d+)?)",  # matches "Version:3.3.0"
    min_version="3.0.0",
    tested_against="3.3.0",
    error_class_missing=PicardNotFoundError,
    is_jar=True,
)

JAVA_SPEC = ToolSpec(
    binary="java",
    version_args=["-version"],  # writes to stderr, not stdout
    version_regex=r'version\s+"(\d+)(?:\.(\d+))?',  # "21.0.10" -> 21; "1.8.0" -> 1
    min_version="11",
    tested_against="21.0.10",
    error_class_missing=JavaNotFoundError,
)

MOSDEPTH_SPEC = ToolSpec(
    binary="mosdepth",
    version_args=["--version"],
    version_regex=r"^mosdepth\s+(\d+\.\d+(?:\.\d+)?)",
    min_version="0.3.6",
    tested_against="0.3.6",
    error_class_missing=MosdepthNotFoundError,
)


def _resolve_picard_jar() -> Path:
    """Find Picard JAR via $PICARD_JAR, conda paths, picard wrapper, or known install dirs.

    Resolution order (B1: conda paths added — most common modern install method):
      1. $PICARD_JAR env var (if set, must point to a real .jar)
      2. ${CONDA_PREFIX}/share/picard-*/picard.jar (bioconda recipe)
         + ${CONDA_PREFIX}/share/picard/picard.jar
         + ${CONDA_PREFIX}/opt/picard-*/picard.jar (alternative recipe layout)
      3. picard wrapper on PATH (script-style installs; we don't parse the script in v0.1)
      4. /usr/share/java/picard.jar, /opt/picard/picard.jar, ~/tools/picard.jar, ~/picard.jar
    """
    env_jar = os.environ.get("PICARD_JAR")
    if env_jar:
        p = Path(env_jar)
        if not p.exists():
            raise PicardNotFoundError(
                what=f"$PICARD_JAR={env_jar}",
                why="PICARD_JAR env var points to nonexistent file",
                fix="Verify $PICARD_JAR or unset and rely on auto-detect",
            )
        return p

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        for pattern in (
            f"{conda_prefix}/share/picard-*/picard.jar",
            f"{conda_prefix}/share/picard/picard.jar",
            f"{conda_prefix}/opt/picard-*/picard.jar",
        ):
            matches = sorted(glob.glob(pattern))
            if matches:
                return Path(matches[-1])  # newest version-suffixed entry

    # picard wrapper on PATH — script-style installs. v0.1 doesn't parse the wrapper
    # to extract the JAR path; if the wrapper is the only thing available, we fall
    # through to the known install paths below.
    if shutil.which("picard"):
        log.debug("`picard` wrapper script on PATH but JAR not yet located; trying known paths")

    for candidate in (
        "/usr/share/java/picard.jar",
        "/opt/picard/picard.jar",
        Path.home() / "tools" / "picard.jar",
        Path.home() / "picard.jar",
    ):
        p = Path(candidate)
        if p.exists():
            return p

    raise PicardNotFoundError(
        what="picard.jar",
        why=(
            "not found via $PICARD_JAR, $CONDA_PREFIX/share/picard-*, picard wrapper, "
            "or known install paths"
        ),
        fix=(
            "Install Picard >= 3.0 via `conda install -c bioconda picard`, or set "
            "$PICARD_JAR=/path/to/picard.jar, or place the JAR at /usr/share/java/picard.jar"
        ),
    )


class ToolWrapper:
    """Subprocess plumbing for one external tool.

    Cached per-spec; binary lookup + version check happen at __init__.
    """

    _version_cache: ClassVar[dict[str, str]] = {}

    def __init__(self, spec: ToolSpec, *, skip_version_check: bool = False):
        self.spec = spec
        self.binary_path = self._resolve_binary(spec)
        if not skip_version_check:
            self._check_version()

    def _resolve_binary(self, spec: ToolSpec) -> Path:
        """Find the binary on PATH (or JAR via _resolve_picard_jar for is_jar specs)."""
        if spec.is_jar:
            return _resolve_picard_jar()
        which = shutil.which(spec.binary)
        if which is None:
            raise spec.error_class_missing(
                what=spec.binary,
                why="binary not found on PATH",
                fix=(
                    f"Install {spec.binary} >= {spec.min_version}; "
                    f"verified against {spec.tested_against} in v2.1 of the spec"
                ),
            )
        return Path(which)

    def _check_version(self) -> None:
        """Probe version, compare against min_version. Raises ToolVersionError on too-old."""
        observed = self.version()
        try:
            obs_v = Version(observed)
            min_v = Version(self.spec.min_version)
        except InvalidVersion as e:
            raise self.spec.error_class_version(
                what=self.spec.binary,
                why=f"unparseable version: observed={observed!r}, min={self.spec.min_version!r}",
                fix=(
                    "This is likely a bug in pileup-aadr's version regex; please file a report. "
                    f"Workaround: install {self.spec.binary} {self.spec.tested_against}"
                ),
            ) from e
        if obs_v < min_v:
            raise self.spec.error_class_version(
                what=self.spec.binary,
                why=f"version too old: observed={observed}, required >= {self.spec.min_version}",
                fix=(
                    f"Upgrade {self.spec.binary} to >= {self.spec.min_version} "
                    f"(verified against {self.spec.tested_against})"
                ),
            )
        log.debug("%s version OK: %s (>= %s)", self.spec.binary, observed, self.spec.min_version)

    def version(self) -> str:
        """Probe the binary for its version string. Cached per-binary.

        Per LLD M15: this intentionally uses subprocess.run directly rather than self.run().
        Two reasons:
          1. Chicken-and-egg — version() is called from __init__ before any tempdir
             exists for stderr-disk discipline.
          2. Bounded output — version probes emit < 1 KB combined stdout+stderr; the
             "no unbounded buffering" discipline rule explicitly carves out small probes.
        """
        cache_key = str(self.binary_path)
        if cache_key in self._version_cache:
            return self._version_cache[cache_key]
        cmd = self._build_invocation([*self.spec.version_args])
        try:
            r = subprocess.run(
                cmd, capture_output=True, text=True, timeout=30, check=False
            )
        except subprocess.TimeoutExpired as e:
            raise self.spec.error_class_version(
                what=self.spec.binary,
                why=f"version probe timed out: {e}",
                fix=f"Verify {self.spec.binary} runs `--version` cleanly",
            ) from e
        combined = (r.stdout or "") + (r.stderr or "")
        m = re.search(self.spec.version_regex, combined, re.MULTILINE)
        if m is None:
            raise self.spec.error_class_version(
                what=self.spec.binary,
                why=f"version regex {self.spec.version_regex!r} matched no line in version output",
                fix=(
                    f"Tool output may have changed. Verify version manually; v2.1 tested against "
                    f"{self.spec.tested_against}"
                ),
            )
        observed = m.group(1)
        self._version_cache[cache_key] = observed
        return observed

    def _build_invocation(
        self, args: list[str], jvm_args: list[str] | None = None
    ) -> list[str]:
        """Construct argv list, prepending `java <jvm_args> -jar <jar>` for is_jar specs.

        For non-jar specs, jvm_args MUST be empty/None (raises ValueError if not).
        """
        if self.spec.is_jar:
            jvm = list(jvm_args or [])
            return ["java", *jvm, "-jar", str(self.binary_path), *args]
        if jvm_args:
            raise ValueError(
                f"jvm_args is only valid for is_jar specs (got for {self.spec.binary})"
            )
        return [str(self.binary_path), *args]

    def run(
        self,
        args: list[str],
        *,
        jvm_args: list[str] | None = None,
        stdin: int | IO | None = None,
        stdout: int | IO | None = None,
        capture_stderr_to: Path,
        check: bool = True,
        timeout: float | None = None,
        extra_env: dict[str, str] | None = None,
    ) -> ToolRunResult:
        """Run the tool; capture stderr to a file.

        For Picard JVM heap (`-Xmx`), pass via jvm_args (NOT in args). The wrapper
        injects them between `java` and `-jar` so the jar receives only its own flags.

        Args:
            args: tool args (without binary name; ToolWrapper adds it)
            jvm_args: JVM flags for is_jar tools (e.g., ["-Xmx3g"]); raises ValueError
                if passed for non-jar specs
            stdin/stdout: subprocess.Popen-compatible
            capture_stderr_to: where to redirect stderr (always written to disk)
            check: if True, raise ToolSubprocessError on non-zero exit
            timeout: kill subprocess if it exceeds N seconds
            extra_env: additional env vars (merged with os.environ)

        Returns:
            ToolRunResult (exit_code, stderr_path, wallclock, peak_rss_mb).
        """
        cmd = self._build_invocation(args, jvm_args=jvm_args)
        env = dict(os.environ)
        if extra_env:
            env.update(extra_env)

        log.debug("Running: %s", " ".join(cmd))

        t0 = time.perf_counter()
        with open(capture_stderr_to, "wb") as stderr_fp:
            try:
                proc = subprocess.run(
                    cmd,
                    stdin=stdin,
                    stdout=stdout,
                    stderr=stderr_fp,
                    timeout=timeout,
                    env=env,
                    check=False,
                )
            except subprocess.TimeoutExpired as e:
                raise ToolSubprocessError(
                    what=self.spec.binary,
                    why=f"timed out after {timeout}s",
                    fix=f"Increase timeout or check for hung subprocess; stderr: {capture_stderr_to}",
                ) from e
        wallclock = time.perf_counter() - t0
        # NOTE: subprocess.run returns CompletedProcess after the child has
        # already exited; /proc/<pid> is gone, so per-child VmHWM read is
        # impossible from this path. RSS measurement is best-effort + reserved
        # for the Popen-based pipe path (Stage 3). v0.1 leaves this None and
        # users who care about RSS use mosdepth (`coverage` subcommand) or an
        # external watcher. See _try_get_child_peak_rss for the read used by
        # the pipe path.
        peak_rss_mb: float | None = None

        result = ToolRunResult(
            exit_code=proc.returncode,
            stdout=proc.stdout.decode() if proc.stdout else None,
            stderr_path=capture_stderr_to,
            stderr_text=None,
            wallclock_seconds=wallclock,
            peak_rss_mb=peak_rss_mb,
        )
        log.info(
            "%s completed: exit=%d wallclock=%.1fs%s",
            self.spec.binary,
            result.exit_code,
            wallclock,
            f" rss={peak_rss_mb:.0f}MB" if peak_rss_mb else "",
        )
        if check and proc.returncode != 0:
            raise ToolSubprocessError(
                what=f"{self.spec.binary} (stderr at {capture_stderr_to})",
                why=f"exit code {proc.returncode}",
                fix=f"Check stderr at {capture_stderr_to} for tool diagnostic",
            )
        return result

    def pipe(
        self,
        downstream: ToolWrapper,
        *,
        upstream_args: list[str],
        downstream_args: list[str],
        upstream_stderr_to: Path,
        downstream_stderr_to: Path,
    ) -> tuple[ToolRunResult, ToolRunResult]:
        """Spawn `self | downstream` as a subprocess pipe.

        Used for Stage 3: `samtools mpileup | pileupCaller`. SIGPIPE handling:
        upstream gets EPIPE if downstream dies first → exits with code 141 (128 + 13).
        Caller checks downstream first (its output matters); upstream 141 is tolerable
        IFF downstream's exit_code is 0.

        Returns:
            (upstream_result, downstream_result)
        """
        upstream_cmd = self._build_invocation(upstream_args)
        downstream_cmd = downstream._build_invocation(downstream_args)

        log.debug("Piping: %s | %s", " ".join(upstream_cmd), " ".join(downstream_cmd))

        t0 = time.perf_counter()
        with (
            open(upstream_stderr_to, "wb") as up_err,
            open(downstream_stderr_to, "wb") as down_err,
        ):
            up_proc = subprocess.Popen(
                upstream_cmd, stdout=subprocess.PIPE, stderr=up_err
            )
            down_proc = subprocess.Popen(
                downstream_cmd, stdin=up_proc.stdout, stderr=down_err
            )
            # Allow upstream to receive SIGPIPE if downstream exits early
            if up_proc.stdout is not None:
                up_proc.stdout.close()

            down_exit = down_proc.wait()
            up_exit = up_proc.wait()
        wallclock = time.perf_counter() - t0

        log.info(
            "%s | %s completed: upstream_exit=%d downstream_exit=%d wallclock=%.1fs",
            self.spec.binary,
            downstream.spec.binary,
            up_exit,
            down_exit,
            wallclock,
        )

        return (
            ToolRunResult(
                exit_code=up_exit,
                stdout=None,
                stderr_path=upstream_stderr_to,
                stderr_text=None,
                wallclock_seconds=wallclock,
                peak_rss_mb=None,
            ),
            ToolRunResult(
                exit_code=down_exit,
                stdout=None,
                stderr_path=downstream_stderr_to,
                stderr_text=None,
                wallclock_seconds=wallclock,
                peak_rss_mb=None,
            ),
        )


def _try_get_child_peak_rss(pid: int) -> float | None:
    """Best-effort peak-RSS for a child process. Returns None if unavailable.

    Linux: read /proc/<pid>/status:VmHWM. macOS: not available per-process via /proc.
    """
    if sys.platform == "linux":
        try:
            with open(f"/proc/{pid}/status") as f:
                for line in f:
                    if line.startswith("VmHWM:"):
                        kb = int(line.split()[1])
                        return kb / 1024.0
        except (FileNotFoundError, ProcessLookupError):
            return None
    return None


def _get_self_rss_mb() -> float:
    """Current RSS of this process (orchestrator), in MB."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Linux ru_maxrss is in KB; macOS in bytes
    return rss / (1024.0 if sys.platform == "linux" else 1024.0 * 1024.0)


__all__ = [
    "JAVA_SPEC",
    "MOSDEPTH_SPEC",
    "PICARD_SPEC",
    "PILEUPCALLER_SPEC",
    "SAMTOOLS_SPEC",
    "ToolRunResult",
    "ToolSpec",
    "ToolWrapper",
]
