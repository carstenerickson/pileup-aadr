"""Reproducible end-to-end perf bench for the `extract` pipeline.

Per LLD §"Performance benchmark (test #20)". Wallclock is hardware-dependent,
so the test reports the elapsed time + asserts a generous ceiling (10 min on
default chr22 fixture). Tighter regression detection is the user's job — pin
your own threshold once you've established a baseline on your hardware.

The benchmark intentionally points at user-supplied bench data via env vars:

    PILEUP_AADR_BENCH_BAM   path to a BAM (chr22-only or whole-genome)
    PILEUP_AADR_BENCH_SNP   path to an AADR .snp slice (any size)
    PILEUP_AADR_BENCH_REF   target FASTA matching the BAM's build

Skipped if any env var is unset or any required external binary is missing.
This avoids shipping multi-GB fixtures in the repo while keeping the bench
runnable for anyone with their own data on disk.

Run via:

    pytest benchmarks/ -v -s -m slow

`-s` so the wallclock readout reaches stderr.
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import pytest

_BAM_ENV = "PILEUP_AADR_BENCH_BAM"
_SNP_ENV = "PILEUP_AADR_BENCH_SNP"
_REF_ENV = "PILEUP_AADR_BENCH_REF"
_BUDGET_SEC = 10 * 60  # 10 min — generous; tighter is hardware-dependent


def _bench_paths_or_skip() -> tuple[Path, Path, Path]:
    missing = [k for k in (_BAM_ENV, _SNP_ENV, _REF_ENV) if not os.environ.get(k)]
    if missing:
        pytest.skip(
            f"Bench env vars unset: {', '.join(missing)} (set to point at "
            "your own data; see benchmarks/bench_extract.py docstring)",
        )
    return tuple(Path(os.environ[k]) for k in (_BAM_ENV, _SNP_ENV, _REF_ENV))  # type: ignore[return-value]


def _has(name: str) -> bool:
    import shutil
    return shutil.which(name) is not None


@pytest.mark.slow
@pytest.mark.skipif(
    not (_has("samtools") and _has("pileupCaller")),
    reason="bench requires samtools + pileupCaller on PATH",
)
def test_extract_perf(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Run extract end-to-end against user-supplied bench data; report wallclock."""
    from pileup_aadr.extract_orch import run_extract
    from pileup_aadr.types import ExtractCliArgs

    bam, snp, ref = _bench_paths_or_skip()
    args = ExtractCliArgs(
        bam=bam,
        aadr_snp=snp,
        ref_fasta=ref,
        output_prefix=tmp_path / "bench_out",
        seed=42,
        # Relaxed gates: bench data may not hit production-coverage levels
        liftover_yield_fail_pct=1.0,
        min_coverage=0,
        warn_coverage=0,
    )
    t0 = time.perf_counter()
    exit_code = run_extract(args)
    elapsed = time.perf_counter() - t0
    assert exit_code == 0

    # Print to stderr so `pytest -s` surfaces the readout regardless of capture.
    with capsys.disabled():
        print(  # noqa: T201 — bench readout is the whole point
            f"\n[bench] extract end-to-end: {elapsed:.1f}s "
            f"(BAM={bam.name}, SNP={snp.name})",
        )

    assert elapsed < _BUDGET_SEC, (
        f"Perf regression: {elapsed:.1f}s > {_BUDGET_SEC}s ceiling. "
        "Set tighter ceilings via your own pytest.ini or by overriding _BUDGET_SEC."
    )
