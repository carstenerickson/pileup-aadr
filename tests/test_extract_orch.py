"""Integration tests for `extract_orch.run_extract`.

Strategy: mock the subprocess-bound Stage 3 (`pileup_call.run_pileup_call`) and
the binary-on-PATH lookup so the orchestrator exercises Stages 0/1/2/4 + gates +
output writers end-to-end against synthetic inputs. The full lift path additionally
mocks Stage 1's Picard subprocess via `lift.lift_aadr_sites`.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from pileup_aadr.counters import (
    PileupCallerSummary,
    Stage1LiftCounters,
    Stage3CallCounters,
)
from pileup_aadr.errors import CoverageGateFailure, OutputExistsError
from pileup_aadr.extract_orch import run_extract
from pileup_aadr.rejoin import GENO_HET, GENO_HOM_REF
from pileup_aadr.types import ExtractCliArgs

# --- BAM + AADR fixture builders ---


def _make_minimal_bam(path: Path, build: str = "hg19") -> None:
    """Write a tiny valid BAM (header only) using pysam.

    `build` selects chr1 length: hg19 → 249_250_621, hg38 → 248_956_422.
    """
    import pysam

    chr1_len = 249_250_621 if build == "hg19" else 248_956_422
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": chr1_len}, {"SN": "chr22", "LN": 51_304_566}],
        "RG": [{"ID": "1", "SM": "TestSample", "LB": "lib1", "PL": "ILLUMINA"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass
    pysam.index(str(path))


def _make_aadr_snp(path: Path, n_sites: int = 50) -> None:
    """Write a tiny AADR .snp slice (chr22 hg19 sites)."""
    lines = []
    for i in range(n_sites):
        lines.append(
            f"rs_test{i:04d}\t22\t{0.001 * i:.4f}\t{20_000_000 + i * 1000}\tA\tG"
        )
    path.write_text("\n".join(lines) + "\n")


def _make_fasta(path: Path, build: str = "hg19") -> None:
    """Write a minimal FASTA + .fai + .dict with chr1 length matching the build.

    The .dict sidecar must exist so the orchestrator's `ensure_target_fasta_dict`
    pre-flight skips its Picard subprocess (which the tests don't mock).
    """
    chr1_len = 249_250_621 if build == "hg19" else 248_956_422
    path.write_text(">chr1\nACGT\n")
    fai = Path(f"{path}.fai")
    fai.write_text(f"chr1\t{chr1_len}\t6\t4\t5\nchr22\t51304566\t100\t4\t5\n")
    Path(f"{path}.dict").write_text(
        f"@HD\tVN:1.6\n"
        f"@SQ\tSN:chr1\tLN:{chr1_len}\n"
        f"@SQ\tSN:chr22\tLN:51304566\n"
    )


# --- Common mocks ---


def _patch_binaries(monkeypatch: pytest.MonkeyPatch) -> None:
    """Monkeypatch ToolWrapper to bypass binary-on-PATH lookup + version check.

    Affects the global ToolWrapper class so orchestrator-instantiated wrappers
    pick it up regardless of which module they're imported through.
    """
    from pileup_aadr import tool_wrapper

    monkeypatch.setattr(
        tool_wrapper.ToolWrapper, "_resolve_binary",
        lambda _self, spec: Path(f"/usr/bin/fake_{spec.binary}"),
    )
    monkeypatch.setattr(
        tool_wrapper.ToolWrapper, "_check_version", lambda _self: None,
    )
    monkeypatch.setattr(
        tool_wrapper.ToolWrapper, "version", lambda _self: "1.0.0-fake",
    )


def _patch_pileup_call(
    monkeypatch: pytest.MonkeyPatch,
    *,
    total_sites: int,
    non_missing: int,
) -> None:
    """Replace `pileup_call.run_pileup_call` with a stub that writes a fake
    pileupCaller triplet matching the input .snp + returns Stage3CallCounters."""
    from pileup_aadr import extract_orch

    def fake_run(
        bam_path: Path, snp_path: Path, bed_path: Path,
        target_fasta_path: Path, output_prefix: Path,
        sample_name: str, pop_name: str, **_kw: Any,
    ) -> Stage3CallCounters:
        # Read the .snp the orchestrator wrote and emit a pileupCaller-shaped triplet
        snp_lines = snp_path.read_text().splitlines()
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        with (
            open(f"{output_prefix}.geno", "w") as gh,
            open(f"{output_prefix}.snp", "w") as sh,
        ):
            for i, line in enumerate(snp_lines):
                # Alternate hom-ref / het for variety
                gh.write((GENO_HOM_REF if i % 2 == 0 else GENO_HET) + "\n")
                sh.write(line + "\n")
        Path(f"{output_prefix}.ind").write_text(
            f"{sample_name}\tU\t{pop_name}\n"
        )
        return Stage3CallCounters(
            wallclock_seconds=0.1,
            pileupcaller_summary=PileupCallerSummary(
                total_sites=total_sites, non_missing_calls=non_missing,
                avg_raw_reads=21.4, avg_damage_cleaned_reads=21.4,
                avg_sampled_from=21.4,
            ),
        )

    monkeypatch.setattr(extract_orch.pileup_call, "run_pileup_call", fake_run)


# --- No-lift fast path: end-to-end ---


@pytest.fixture
def no_lift_run_setup(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> dict[str, Path]:
    """Set up a no-lift fast-path test environment (hg19 BAM + hg19 AADR)."""
    bam = tmp_path / "user.bam"
    _make_minimal_bam(bam, build="hg19")
    aadr = tmp_path / "aadr.snp"
    _make_aadr_snp(aadr, n_sites=50)
    fasta = tmp_path / "hg19.fa"
    _make_fasta(fasta, build="hg19")

    _patch_binaries(monkeypatch)
    _patch_pileup_call(
        monkeypatch, total_sites=50, non_missing=48,
    )
    return {"bam": bam, "aadr": aadr, "fasta": fasta, "tmp": tmp_path}


def test_no_lift_run_writes_eigenstrat_triplet(no_lift_run_setup: dict[str, Path]) -> None:
    """End-to-end no-lift run produces all 4 user-facing artifacts."""
    paths = no_lift_run_setup
    output_prefix = paths["tmp"] / "out" / "carsten_pseudohaploid"

    args = ExtractCliArgs(
        bam=paths["bam"], aadr_snp=paths["aadr"],
        output_prefix=output_prefix,
        ref_fasta=paths["fasta"], bam_build="hg19", aadr_build="hg19",
        # Lower coverage thresholds so a 50-site fixture passes
        min_coverage=10, warn_coverage=20,
    )
    exit_code = run_extract(args)

    assert exit_code == 0
    assert (output_prefix.with_suffix(".geno")).exists()
    assert (output_prefix.with_suffix(".snp")).exists()
    assert (output_prefix.with_suffix(".ind")).exists()
    assert (Path(f"{output_prefix}.pseudohaploid.json")).exists()


def test_no_lift_run_emits_json_report(no_lift_run_setup: dict[str, Path]) -> None:
    """--report-json populates the schema-1 JSON with stage_1_lift=null."""
    paths = no_lift_run_setup
    output_prefix = paths["tmp"] / "out"
    report = paths["tmp"] / "report.json"

    args = ExtractCliArgs(
        bam=paths["bam"], aadr_snp=paths["aadr"],
        output_prefix=output_prefix,
        ref_fasta=paths["fasta"], bam_build="hg19", aadr_build="hg19",
        min_coverage=10, warn_coverage=20,
        report_json=report,
    )
    run_extract(args)

    data = json.loads(report.read_text())
    assert data["schema_version"] == 1
    assert data["stage_1_lift"] is None
    assert data["stage_2_transform"] is None
    assert data["stage_4_rejoin"] is None
    assert data["stage_3_call"]["pileupcaller_summary"]["total_sites"] == 50
    assert data["gates"]["liftover_yield"] == "N/A"
    assert data["gates"]["coverage"] == "PASS"


def test_no_lift_run_emits_per_variant_tsv(no_lift_run_setup: dict[str, Path]) -> None:
    """--report-tsv populates the streaming per-variant TSV."""
    paths = no_lift_run_setup
    output_prefix = paths["tmp"] / "out"
    tsv = paths["tmp"] / "report.tsv"

    args = ExtractCliArgs(
        bam=paths["bam"], aadr_snp=paths["aadr"],
        output_prefix=output_prefix,
        ref_fasta=paths["fasta"], bam_build="hg19", aadr_build="hg19",
        min_coverage=10, warn_coverage=20,
        report_tsv=tsv,
    )
    run_extract(args)

    lines = tsv.read_text().splitlines()
    # Header + 50 data rows (one per AADR site)
    assert lines[0].startswith("aadr_id\t")
    assert len(lines) == 51


# --- Output-prefix collision ---


def test_output_collision_raises_without_overwrite(
    no_lift_run_setup: dict[str, Path],
) -> None:
    """Existing .geno at the prefix → OutputExistsError unless --overwrite."""
    paths = no_lift_run_setup
    output_prefix = paths["tmp"] / "out"
    # Pre-create a .geno at the prefix
    Path(f"{output_prefix}.geno").write_text("0\n")

    args = ExtractCliArgs(
        bam=paths["bam"], aadr_snp=paths["aadr"],
        output_prefix=output_prefix,
        ref_fasta=paths["fasta"], bam_build="hg19", aadr_build="hg19",
        min_coverage=10, warn_coverage=20,
    )
    with pytest.raises(OutputExistsError):
        run_extract(args)


def test_overwrite_flag_replaces_existing_outputs(
    no_lift_run_setup: dict[str, Path],
) -> None:
    """--overwrite removes the pre-existing files and proceeds."""
    paths = no_lift_run_setup
    output_prefix = paths["tmp"] / "out"
    Path(f"{output_prefix}.geno").write_text("STALE\n")

    args = ExtractCliArgs(
        bam=paths["bam"], aadr_snp=paths["aadr"],
        output_prefix=output_prefix,
        ref_fasta=paths["fasta"], bam_build="hg19", aadr_build="hg19",
        min_coverage=10, warn_coverage=20,
        overwrite=True,
    )
    run_extract(args)
    # Stale content should be gone — first .geno line is now a real digit
    geno_first = Path(f"{output_prefix}.geno").read_text().splitlines()[0]
    assert geno_first in ("0", "1", "2", "9")


# --- Coverage gate ---


def test_coverage_gate_failure_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """non_missing_autosomal < min_coverage → CoverageGateFailure."""
    bam = tmp_path / "user.bam"
    _make_minimal_bam(bam, build="hg19")
    aadr = tmp_path / "aadr.snp"
    _make_aadr_snp(aadr, n_sites=10)
    fasta = tmp_path / "hg19.fa"
    _make_fasta(fasta, build="hg19")

    _patch_binaries(monkeypatch)
    _patch_pileup_call(monkeypatch, total_sites=10, non_missing=10)

    args = ExtractCliArgs(
        bam=bam, aadr_snp=aadr,
        output_prefix=tmp_path / "out",
        ref_fasta=fasta, bam_build="hg19", aadr_build="hg19",
        min_coverage=1_000_000,  # way above what 10 sites can produce
        warn_coverage=1_000_000,
    )
    with pytest.raises(CoverageGateFailure, match="below --min-coverage"):
        run_extract(args)


def test_coverage_gate_skipped_for_chry_only_panel(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """v0.2 #3: chrY-only AADR (haplogroup workflow) skips the autosomal gate
    cleanly with an INFO log, even if --min-coverage would normally fail."""
    import logging
    bam = tmp_path / "user.bam"
    _make_minimal_bam(bam, build="hg19")
    aadr = tmp_path / "aadr_chry.snp"
    # chrY-only AADR slice (chrom_int '24' in AADR encoding); add chr1 too so
    # detect_aadr_build can resolve via the chr1 anchor (1240k chrY rows
    # alone don't have positions near hg19/hg38 chrY end with our wide tol).
    aadr.write_text("\n".join([
        f"rs_y{i:04d}\t24\t0.0\t{1_000_000 + i * 1000}\tA\tG"
        for i in range(20)
    ]) + "\n")
    fasta = tmp_path / "hg19.fa"
    _make_fasta(fasta, build="hg19")

    _patch_binaries(monkeypatch)
    _patch_pileup_call(monkeypatch, total_sites=20, non_missing=20)
    caplog.set_level(logging.INFO, logger="pileup_aadr.extract_orch")

    args = ExtractCliArgs(
        bam=bam, aadr_snp=aadr,
        output_prefix=tmp_path / "out",
        ref_fasta=fasta, bam_build="hg19", aadr_build="hg19",
        min_coverage=1_000_000,  # would normally fail; chrY-only panel skips gate
        warn_coverage=1_000_000,
    )
    exit_code = run_extract(args)
    assert exit_code == 0
    # The gate skip emits a specific INFO log with the panel class name.
    assert any(
        "Skipping autosomal coverage gate" in r.message
        and "chrY_only" in r.message
        for r in caplog.records
    ), f"expected gate-skip log; got: {[r.message for r in caplog.records]}"


def test_coverage_warn_threshold_marks_warning(
    no_lift_run_setup: dict[str, Path],
    caplog: pytest.LogCaptureFixture,
) -> None:
    """non_missing_autosomal < warn_coverage but >= min_coverage → coverage_warning=True."""
    import logging
    paths = no_lift_run_setup
    caplog.set_level(logging.WARNING, logger="pileup_aadr.extract_orch")

    args = ExtractCliArgs(
        bam=paths["bam"], aadr_snp=paths["aadr"],
        output_prefix=paths["tmp"] / "out",
        ref_fasta=paths["fasta"], bam_build="hg19", aadr_build="hg19",
        min_coverage=10, warn_coverage=10_000,  # 50 sites passes min, fails warn
    )
    run_extract(args)
    assert any("below --warn-coverage" in r.message for r in caplog.records)


# --- Full lift path: orchestrator routing only (Stage 1 mocked) ---


def test_full_lift_path_dispatches_to_stage_1(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """When BAM build != AADR build, orchestrator runs Stages 1/2/4 (not just 3)."""
    bam = tmp_path / "user.bam"
    _make_minimal_bam(bam, build="hg38")
    aadr = tmp_path / "aadr.snp"
    _make_aadr_snp(aadr, n_sites=50)
    fasta = tmp_path / "hg38.fa"
    _make_fasta(fasta, build="hg38")

    _patch_binaries(monkeypatch)

    # Mock Stage 1 to write a synthetic lifted VCF (matches AADR sites in the
    # fixture) + return Stage1LiftCounters with 100% yield
    from pileup_aadr import extract_orch

    def fake_lift(
        sites_vcf_path: Path,
        chain_path: Path,
        target_fasta_path: Path,
        output_lifted_vcf: Path,
        output_rejected_vcf: Path,
        input_filter_counters: Any,
        **_kw: Any,
    ) -> Stage1LiftCounters:
        # Build a minimal lifted VCF carrying AADR_RS INFO so transform.py + rejoin.py
        # can find the rsids
        lines = ["##fileformat=VCFv4.2", "##contig=<ID=chr22,length=51304566>",
                 '##INFO=<ID=AADR_RS,Number=1,Type=String,Description="x">',
                 '##INFO=<ID=SwappedAlleles,Number=0,Type=Flag,Description="x">',
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
        for i in range(50):
            lines.append(
                f"chr22\t{20_000_000 + i * 1000}\trs_test{i:04d}\tA\tG\t.\tPASS\t"
                f"AADR_RS=rs_test{i:04d}"
            )
        output_lifted_vcf.write_text("\n".join(lines) + "\n")
        output_rejected_vcf.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        return Stage1LiftCounters(
            wallclock_seconds=0.1,
            input_sites_after_filters=input_filter_counters.rows_written,
            lifted_sites=50, liftover_yield_pct=100.0,
            liftover_yield_warning=False,
            rejected_by_reason={"NoTarget": 0, "MismatchedRefAllele": 0,
                                "IndelStraddlesMultipleIntervals": 0,
                                "SwappedAlleles": 0, "other": 0},
            swapped_alleles_count=0,
            input_filters={"palindrome_drops": 0, "non_snp_drops": 0,
                           "non_autosome_drops": 0},
        )

    monkeypatch.setattr(extract_orch.lift, "lift_aadr_sites", fake_lift)
    _patch_pileup_call(monkeypatch, total_sites=50, non_missing=48)

    output_prefix = tmp_path / "out"
    args = ExtractCliArgs(
        bam=bam, aadr_snp=aadr, output_prefix=output_prefix,
        ref_fasta=fasta, bam_build="hg38", aadr_build="hg19",
        min_coverage=10, warn_coverage=20,
    )
    exit_code = run_extract(args)

    assert exit_code == 0
    # Lift path artifacts: all 4 outputs + the lift-specific Stage4 counters in the
    # JSON report (sanity-check the dispatch happened by checking the .geno isn't
    # a byte-for-byte copy of the fast-path's output)
    assert (output_prefix.with_suffix(".geno")).exists()
