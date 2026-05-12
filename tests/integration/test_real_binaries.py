"""End-to-end tests against real samtools / pileupCaller / picard / mosdepth.

These tests exist to catch the class of bugs the mocked unit suite can't see —
exactly the three smoke-test bugs caught on ancestrytracke-f (Picard --version
regex, CompletedProcess.pid AttributeError, scientific-notation parser drift).

Each test gates on the corresponding `requires_*` marker and skips silently
when the binary isn't installed. CI's bio-tools job installs all four via
bioconda + sets PICARD_JAR; local `pip install -e '.[dev]'` skips them.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pileup_aadr.extract_orch import run_extract
from pileup_aadr.tool_wrapper import (
    MOSDEPTH_SPEC,
    PICARD_SPEC,
    PILEUPCALLER_SPEC,
    SAMTOOLS_SPEC,
    ToolWrapper,
)
from pileup_aadr.types import CoverageCliArgs, ExtractCliArgs

from .conftest import (
    requires_mosdepth,
    requires_picard,
    requires_pileupcaller,
    requires_samtools,
)

# --- Version probes (catch the Picard --version regression) ---


@requires_samtools
def test_samtools_version_probe() -> None:
    """ToolWrapper(SAMTOOLS_SPEC).version() returns a parseable version string."""
    v = ToolWrapper(SAMTOOLS_SPEC).version()
    assert v
    assert "." in v


@requires_pileupcaller
def test_pileupcaller_version_probe() -> None:
    v = ToolWrapper(PILEUPCALLER_SPEC).version()
    assert v
    assert "." in v


@requires_picard
def test_picard_version_probe() -> None:
    """Regression: PICARD_SPEC.version_args=['LiftoverVcf', '--version'] (the
    no-subcommand --version prints help, not the version string)."""
    v = ToolWrapper(PICARD_SPEC).version()
    assert v
    assert "." in v


@requires_mosdepth
def test_mosdepth_version_probe() -> None:
    v = ToolWrapper(MOSDEPTH_SPEC).version()
    assert v
    assert "." in v


# --- End-to-end extract (catches CompletedProcess.pid + parser drift) ---


@requires_samtools
@requires_pileupcaller
@requires_picard
@pytest.mark.slow
def test_extract_end_to_end(
    hg38_chr22_fasta: Path,
    hg38_chr22_bam: Path,
    aadr_chr22_slice: Path,
    tmp_path: Path,
) -> None:
    """Full pipeline against real binaries: AADR slice (hg19) + chr22 BAM (hg38).

    Empty BAM → 0 non-missing calls but a parseable pileupCaller summary.
    Validates: Stage 1 lift exits 0 + parses; Stage 2 transform; Stage 3 pipe
    + scientific-notation-parser; Stage 4 rejoin; coverage gate (relaxed to
    accept zero); JSON report writes.
    """
    out_prefix = tmp_path / "out" / "smoke"
    report = tmp_path / "smoke.report.json"

    args = ExtractCliArgs(
        bam=hg38_chr22_bam,
        aadr_snp=aadr_chr22_slice,
        output_prefix=out_prefix,
        ref_fasta=hg38_chr22_fasta,
        bam_build="hg38",
        aadr_build="hg19",
        picard_mem="2g",
        report_json=report,
        # Relaxed gates: 50-site slice with empty BAM produces ~0 calls
        liftover_yield_fail_pct=1.0,
        liftover_yield_warn_pct=50.0,
        min_coverage=0,
        warn_coverage=0,
    )
    exit_code = run_extract(args)
    assert exit_code == 0

    # All four output artifacts present
    for ext in (".geno", ".snp", ".ind", ".pseudohaploid.json"):
        assert Path(f"{out_prefix}{ext}").exists(), f"missing output {ext}"

    # JSON report has the expected schema-1 shape
    data = json.loads(report.read_text())
    assert data["schema_version"] == 1
    assert data["stage_1_lift"] is not None
    assert data["stage_3_call"]["pileupcaller_summary"]["total_sites"] >= 0
    # Load-bearing invariant: Stage1 swap count == Stage4 swap count
    assert (
        data["stage_1_lift"]["swapped_alleles_count"]
        == data["stage_4_rejoin"]["ref_alt_swap_count"]
    )


# --- Coverage subcommand (mosdepth real run) ---


@requires_mosdepth
@requires_samtools
def test_coverage_against_real_mosdepth(
    hg38_chr22_bam: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """`coverage` subcommand emits the 9-col TSV against a real mosdepth run."""
    from pileup_aadr.coverage_impl import run_coverage

    args = CoverageCliArgs(bam=hg38_chr22_bam)
    exit_code = run_coverage(args)
    assert exit_code == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0] == (
        "chrom\tlength\tbases\tmean_coverage\tmedian_coverage"
        "\tfraction_at_>=1x\tfraction_at_>=5x"
        "\tfraction_at_>=10x\tfraction_at_>=30x"
    )
    # Empty BAM → mosdepth emits only the `total` rollup row (no per-chrom rows
    # without reads). The CI bio-tools job verifies the schema; real BAMs would
    # add per-chrom rows for any chrom with reads.
    assert any(line.startswith("total\t") for line in out)
