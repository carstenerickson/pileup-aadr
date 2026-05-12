"""LLD test #19: pileup-aadr → pgen-samplebind merge → AT2 extract_f2 stays
stable to within `max_dev < 1e-9` of a frozen f2 baseline.

Two layers, both gated on AT2 (R + admixtools) + pgen-samplebind:

- **Layer A** (always runs when toolchain present): loads a pre-extracted
  pileup-aadr triplet from `fixtures/lld19/pileup_aadr_output/`, merges with
  the AADR chr22 5-pop reference subset, runs AT2 `extract_f2`, asserts the
  result matches `fixtures/lld19/baseline/f2_baseline.rds`. Validates the
  (samplebind merge → AT2 extract_f2) composition stays stable. ~10 sec.

- **Layer B** (skipped without `PILEUP_AADR_LLD19_BAM`): full chain — runs
  `pileup-aadr extract` against the env-var BAM, then merges + extract_f2
  + diff baseline. Validates the WHOLE chain including pileup-aadr drift.
  ~1-2 min on a chr22-subset BAM.

The baseline was generated on the `ancestrytracke-f` reference instance
with pileup-aadr v0.1.2 + pgen-samplebind v0.1.0 + AT2 production/v1.0.
See `fixtures/lld19/README.md` and `scripts/regen_lld19_baseline.sh` for
regeneration procedure.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

LLD19_DIR = Path(__file__).parent / "fixtures" / "lld19"
AADR_REF_PREF = LLD19_DIR / "aadr_ref" / "aadr_chr22_5pops"
PRE_EXTRACTED_PREF = LLD19_DIR / "pileup_aadr_output" / "extract_out"
BASELINE_RDS = LLD19_DIR / "baseline" / "f2_baseline.rds"

# Tolerance per LLD §"pgen-samplebind composition" #19.
MAX_DEV_THRESHOLD: float = 1e-9


def _has_admixtools() -> bool:
    """True iff R + admixtools are importable. Validates the AT2 fork is installed."""
    if shutil.which("Rscript") is None:
        return False
    proc = subprocess.run(
        ["Rscript", "-e",
         "if (requireNamespace('admixtools', quietly=TRUE)) cat('ok')"],
        capture_output=True, text=True, check=False,
    )
    return proc.returncode == 0 and "ok" in proc.stdout


def _has_pgen_samplebind() -> bool:
    return shutil.which("pgen-samplebind") is not None


requires_admixtools = pytest.mark.skipif(
    not _has_admixtools(),
    reason="AT2 not installed (Rscript -e 'library(admixtools)' must succeed)",
)
requires_pgen_samplebind = pytest.mark.skipif(
    not _has_pgen_samplebind(),
    reason="pgen-samplebind binary not on PATH",
)


def _run_pgen_samplebind_merge(
    aadr_pref: Path, target_pref: Path, out_pref: Path,
) -> None:
    """Run pgen-samplebind merge: AADR ref + target → merged PFILE at out_pref."""
    proc = subprocess.run(
        [
            "pgen-samplebind", "merge",
            str(aadr_pref),
            "--target", str(target_pref),
            "-o", str(out_pref),
            "--quiet",
        ],
        capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 0, (
        f"pgen-samplebind merge failed: {proc.stderr}"
    )


def _run_at2_extract_f2_and_diff(
    merged_pref: Path,
    work_dir: Path,
    baseline_rds: Path,
) -> float:
    """Run AT2 extract_f2 on the merged PFILE; return the max abs deviation
    from the baseline f2_blocks array.

    The R script writes the per-pair f2 directory under work_dir/f2_dir,
    loads it via `f2_from_precomp`, loads the baseline .rds, computes the
    elementwise abs-diff max, and prints it to stdout for the Python
    caller to parse.
    """
    r_script = (
        "suppressPackageStartupMessages(library(admixtools))\n"
        f"extract_f2('{merged_pref}', outdir='{work_dir / 'f2_dir'}', "
        "maxmiss=0, overwrite=TRUE, verbose=FALSE, auto_only=FALSE)\n"
        f"test <- f2_from_precomp('{work_dir / 'f2_dir'}', verbose=FALSE)\n"
        f"baseline <- readRDS('{baseline_rds}')\n"
        # Both are 3D arrays [pop x pop x block]. Order MUST match — the test
        # baseline was generated from the same ref pops + sample.
        "stopifnot(identical(dim(test), dim(baseline)))\n"
        "stopifnot(identical(dimnames(test), dimnames(baseline)))\n"
        "max_dev <- max(abs(test - baseline), na.rm=TRUE)\n"
        "cat(sprintf('MAX_DEV=%.20e\\n', max_dev))\n"
    )
    proc = subprocess.run(
        ["Rscript", "-e", r_script],
        capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 0, (
        f"AT2 extract_f2 failed:\nSTDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}"
    )
    # Parse MAX_DEV=... from stdout
    for line in proc.stdout.splitlines():
        if line.startswith("MAX_DEV="):
            return float(line.split("=", 1)[1])
    raise AssertionError(
        f"R script returned no MAX_DEV line. Output:\n{proc.stdout}",
    )


# --- Layer A: composition stability (no BAM needed) ---


@requires_admixtools
@requires_pgen_samplebind
def test_lld19_composition_matches_baseline(tmp_path: Path) -> None:
    """Layer A: pre-extracted pileup-aadr triplet → pgen-samplebind merge →
    AT2 extract_f2 must match the frozen f2 baseline within 1e-9.

    Failure here means either pgen-samplebind or AT2 changed behavior
    on identical input — pinning either via lockfile or version-bumping
    the baseline is the right response."""
    merged_pref = tmp_path / "merged"
    _run_pgen_samplebind_merge(
        aadr_pref=AADR_REF_PREF,
        target_pref=PRE_EXTRACTED_PREF,
        out_pref=merged_pref,
    )
    max_dev = _run_at2_extract_f2_and_diff(
        merged_pref=merged_pref,
        work_dir=tmp_path,
        baseline_rds=BASELINE_RDS,
    )
    assert max_dev < MAX_DEV_THRESHOLD, (
        f"f2 drifted from baseline: max_dev={max_dev:.3e} >= {MAX_DEV_THRESHOLD:.0e}. "
        "Either pgen-samplebind/AT2 changed behavior on identical input, or "
        "the fixtures need regeneration via scripts/regen_lld19_baseline.sh."
    )


# --- Layer B: full chain (gated on a real chr22 BAM via env var) ---


_BAM_ENV = "PILEUP_AADR_LLD19_BAM"
_REF_ENV = "PILEUP_AADR_LLD19_REF"


@requires_admixtools
@requires_pgen_samplebind
@pytest.mark.skipif(
    not (os.environ.get(_BAM_ENV) and os.environ.get(_REF_ENV)),
    reason=(
        f"Bench env vars unset: set {_BAM_ENV} (chr22-coverage BAM, hg38) and "
        f"{_REF_ENV} (matching FASTA) to run the full-chain layer"
    ),
)
@pytest.mark.slow
def test_lld19_full_chain_matches_baseline(tmp_path: Path) -> None:
    """Layer B: full chain (BAM → pileup-aadr extract → samplebind merge →
    AT2 extract_f2) must match the frozen f2 baseline within 1e-9.

    This is the canonical LLD #19 invariant. Failure here means pileup-aadr
    drifted (a Stage 1/2/3/4 change moved a genotype call) OR pgen-samplebind/
    AT2 changed (caught by layer A independently)."""
    from pileup_aadr.extract_orch import run_extract
    from pileup_aadr.types import ExtractCliArgs

    bam = Path(os.environ[_BAM_ENV])
    ref = Path(os.environ[_REF_ENV])
    aadr_snp = AADR_REF_PREF.with_suffix(".snp")

    extract_pref = tmp_path / "extract_out"
    args = ExtractCliArgs(
        bam=bam,
        aadr_snp=aadr_snp,
        output_prefix=extract_pref,
        ref_fasta=ref,
        bam_build="hg38",
        aadr_build="hg19",
        picard_mem="8g",
        sample_name="GFX0442453",
        pop_name="GFX0442453_test",
        # Gates relaxed enough for a chr22-only run
        min_coverage=100,
        warn_coverage=1000,
        liftover_yield_fail_pct=50.0,
    )
    exit_code = run_extract(args)
    assert exit_code == 0

    merged_pref = tmp_path / "merged"
    _run_pgen_samplebind_merge(
        aadr_pref=AADR_REF_PREF,
        target_pref=extract_pref,
        out_pref=merged_pref,
    )
    max_dev = _run_at2_extract_f2_and_diff(
        merged_pref=merged_pref,
        work_dir=tmp_path,
        baseline_rds=BASELINE_RDS,
    )
    assert max_dev < MAX_DEV_THRESHOLD, (
        f"Full-chain f2 drifted from baseline: max_dev={max_dev:.3e}. "
        "If layer A passes but this fails, the drift is in pileup-aadr; "
        "a Stage 1-4 change moved one or more genotype calls."
    )
