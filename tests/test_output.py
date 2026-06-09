"""Tests for output.py — sidecar JSON, JSON report, per-variant TSV, stdout summary."""
from __future__ import annotations

import io
import json
from pathlib import Path
from unittest.mock import patch

from pileup_aadr.counters import (
    CoverageCounters,
    ExtractCounters,
    PileupCallerSummary,
    Stage1InputFilters,
    Stage1LiftCounters,
    Stage2TransformCounters,
    Stage3CallCounters,
    Stage4RejoinCounters,
)
from pileup_aadr.output import (
    JSON_REPORT_SCHEMA_VERSION,
    PSEUDOHAPLOID_SIDECAR_SCHEMA_VERSION,
    write_json_report,
    write_per_variant_tsv,
    write_pseudohaploid_sidecar,
    write_stdout_summary,
)

# --- Fixture builders ---


def _full_lift_counters() -> ExtractCounters:
    """A canonical ExtractCounters for the lift path (all stages populated)."""
    return ExtractCounters(
        stage_1_lift=Stage1LiftCounters(
            wallclock_seconds=12.3,
            input_sites_after_filters=1_100_000,
            lifted_sites=1_098_500,
            liftover_yield_pct=99.8636,
            liftover_yield_warning=False,
            rejected_by_reason={
                "NoTarget": 1500, "MismatchedRefAllele": 0,
                "IndelStraddlesMultipleIntervals": 0,
                "SwappedAlleles": 0, "other": 0,
            },
            swapped_alleles_count=16,
            input_filters=Stage1InputFilters(
                palindrome_drops=95_000,
                non_snp_drops=0,
                non_autosome_drops=0,
                rows_written=1_100_000,
            ),
        ),
        stage_2_transform=Stage2TransformCounters(
            wallclock_seconds=4.2, alt_contig_drops=120, output_sites=1_098_380,
        ),
        stage_3_call=Stage3CallCounters(
            wallclock_seconds=1832.5,
            pileupcaller_summary=PileupCallerSummary(
                total_sites=1_098_380, non_missing_calls=1_087_204,
                avg_raw_reads=21.4, avg_damage_cleaned_reads=21.4,
                avg_sampled_from=21.4,
            ),
        ),
        stage_4_rejoin=Stage4RejoinCounters(
            wallclock_seconds=8.7,
            rsid_matched=1_087_204, ref_alt_swap_count=16,
            allele_mismatch_drops=0, output_variants=1_087_188,
        ),
        coverage=CoverageCounters(
            non_missing_autosomal_calls=1_050_000,
            coverage_fraction=0.9545,
            coverage_warning=False,
            per_chrom_call_count={f"chr{i}": 50_000 for i in range(1, 23)},
        ),
        gates={"liftover_yield": "PASS", "coverage": "PASS"},
        wallclock_total_seconds=1858.2,
    )


def _no_lift_counters() -> ExtractCounters:
    """ExtractCounters for the no-lift fast path (Stages 1/2/4 are None)."""
    return ExtractCounters(
        stage_1_lift=None, stage_2_transform=None,
        stage_3_call=Stage3CallCounters(
            wallclock_seconds=1500.0,
            pileupcaller_summary=PileupCallerSummary(
                total_sites=1_200_000, non_missing_calls=1_180_000,
                avg_raw_reads=21.4, avg_damage_cleaned_reads=21.4,
                avg_sampled_from=21.4,
            ),
        ),
        stage_4_rejoin=None,
        coverage=CoverageCounters(
            non_missing_autosomal_calls=1_150_000,
            coverage_fraction=0.96,
            coverage_warning=False,
            per_chrom_call_count={f"chr{i}": 52_000 for i in range(1, 23)},
        ),
        gates={"liftover_yield": "N/A", "coverage": "PASS"},
        wallclock_total_seconds=1502.0,
    )


# --- write_pseudohaploid_sidecar ---


def test_sidecar_writes_dict_as_json(tmp_path: Path) -> None:
    """Round-trip sidecar dict → JSON → loaded dict equals original."""
    sidecar = {
        "schema_version": 1,
        "samples": {
            "Carsten": {
                "pseudohaploid": 1, "het_count": 100,
                "non_missing_autosomal_count": 1_000_000,
                "het_rate": 0.0001,
                "source": "pileup-aadr-extract",
                "calling_mode": "randomHaploid",
                "note": "--randomHaploid output is pseudohaploid by construction "
                "(one random read per site → 0% het)",
            },
        },
    }
    path = tmp_path / "out.pseudohaploid.json"
    write_pseudohaploid_sidecar(path, sidecar)
    loaded = json.loads(path.read_text())
    assert loaded == sidecar


def test_sidecar_injects_schema_version_when_missing(tmp_path: Path) -> None:
    """If caller omits schema_version, the writer injects the canonical value."""
    sidecar = {"samples": {"S": {"pseudohaploid": 1}}}
    path = tmp_path / "out.json"
    write_pseudohaploid_sidecar(path, sidecar)
    loaded = json.loads(path.read_text())
    assert loaded["schema_version"] == PSEUDOHAPLOID_SIDECAR_SCHEMA_VERSION


def test_sidecar_creates_parent_dir(tmp_path: Path) -> None:
    """Parent dir is created on demand."""
    path = tmp_path / "deep" / "nested" / "out.json"
    write_pseudohaploid_sidecar(path, {"samples": {}})
    assert path.exists()


# --- write_json_report ---


def test_json_report_includes_schema_version_and_tool_block(tmp_path: Path) -> None:
    """Report has schema_version + tool block + ExtractCounters fields at top level."""
    counters = _full_lift_counters()
    path = tmp_path / "report.json"
    write_json_report(
        path, counters,
        config={"threads": 4},
        tool_versions={"samtools": "1.23.1", "pileupCaller": "1.6.0.0"},
        input_meta={"bam_path": "/data/in.bam", "bam_build": "hg38"},
        output_meta={"prefix": str(tmp_path / "out"), "geno_bytes": 100},
    )
    data = json.loads(path.read_text())

    assert data["schema_version"] == JSON_REPORT_SCHEMA_VERSION
    assert data["tool"]["name"] == "pileup-aadr"
    assert data["tool"]["tool_versions"]["samtools"] == "1.23.1"
    assert data["input"]["bam_path"] == "/data/in.bam"
    assert data["output"]["geno_bytes"] == 100
    assert data["config"] == {"threads": 4}
    # ExtractCounters fields land at the top level (not nested under "counters")
    assert data["stage_1_lift"]["lifted_sites"] == 1_098_500
    assert data["stage_3_call"]["pileupcaller_summary"]["total_sites"] == 1_098_380
    assert data["coverage"]["non_missing_autosomal_calls"] == 1_050_000
    assert data["gates"] == {"liftover_yield": "PASS", "coverage": "PASS"}


def test_json_report_no_lift_fast_path_has_null_stages(tmp_path: Path) -> None:
    """No-lift path: stage_1_lift / stage_2_transform / stage_4_rejoin serialize as null."""
    counters = _no_lift_counters()
    path = tmp_path / "report.json"
    write_json_report(
        path, counters, config={}, tool_versions={},
        input_meta={}, output_meta={},
    )
    data = json.loads(path.read_text())
    assert data["stage_1_lift"] is None
    assert data["stage_2_transform"] is None
    assert data["stage_4_rejoin"] is None
    # Stage 3 is always populated
    assert data["stage_3_call"]["pileupcaller_summary"]["non_missing_calls"] == 1_180_000


# --- write_per_variant_tsv ---


def test_per_variant_tsv_writes_header_and_rows(tmp_path: Path) -> None:
    """TSV starts with the 6-col header; row count = input row count."""
    rows = [
        {"aadr_id": "rs1", "chrom_hg19": "1", "pos_hg19": 1000,
         "ref_hg19": "A", "alt_hg19": "G", "action": "passthrough"},
        {"aadr_id": "rs2", "chrom_hg19": "1", "pos_hg19": 2000,
         "ref_hg19": "C", "alt_hg19": "T", "action": "swap"},
    ]
    path = tmp_path / "report.tsv"
    n = write_per_variant_tsv(path, iter(rows))
    assert n == 2

    lines = path.read_text().splitlines()
    assert lines[0] == "aadr_id\tchrom_hg19\tpos_hg19\tref_hg19\talt_hg19\taction"
    assert lines[1] == "rs1\t1\t1000\tA\tG\tpassthrough"
    assert lines[2] == "rs2\t1\t2000\tC\tT\tswap"


def test_per_variant_tsv_streams_constant_memory(tmp_path: Path) -> None:
    """Iterator (not list) is accepted — large reports don't materialize in memory."""
    def gen():
        for i in range(100):
            yield {
                "aadr_id": f"rs{i}", "chrom_hg19": "1", "pos_hg19": i,
                "ref_hg19": "A", "alt_hg19": "G", "action": "passthrough",
            }
    path = tmp_path / "report.tsv"
    n = write_per_variant_tsv(path, gen())
    assert n == 100


# --- write_stdout_summary ---


def test_stdout_summary_lift_path_renders_all_stages(tmp_path: Path) -> None:
    """Lift-path summary mentions all 4 stages + coverage gate."""
    counters = _full_lift_counters()
    buf = io.StringIO()
    with patch("sys.stdout", buf):
        write_stdout_summary(
            counters, bam_path=Path("/data/in.bam"), bam_format="BAM",
            bam_build="hg38", bam_coverage=21.4,
            aadr_path=Path("/data/aadr.snp"), aadr_total=1_200_000,
            ref_fasta=Path("/data/hg38.fa"), chain_path=Path("/data/chain.gz"),
            output_prefix=tmp_path / "out",
            output_bytes={"geno": 100, "snp": 200, "ind": 50, "pseudohaploid_json": 500},
        )
    text = buf.getvalue()
    assert "Stage 1" in text
    assert "Stage 2" in text
    assert "Stage 3" in text
    assert "Stage 4" in text
    assert "Coverage report" in text
    assert "Done in" in text
    # Per-chrom 6-row grid for autosomes
    assert "chr1:" in text
    assert "chr22:" in text


def test_stdout_summary_no_lift_skips_stages_124(tmp_path: Path) -> None:
    """No-lift summary has no Stage 1/2/4 sections; Stage 3 + coverage stay."""
    counters = _no_lift_counters()
    buf = io.StringIO()
    with patch("sys.stdout", buf):
        write_stdout_summary(
            counters, bam_path=Path("/data/in.bam"), bam_format="BAM",
            bam_build="hg19", bam_coverage=None,
            aadr_path=Path("/data/aadr.snp"), aadr_total=1_200_000,
            ref_fasta=Path("/data/hg19.fa"), chain_path=Path("(bundled)"),
            output_prefix=tmp_path / "out",
            output_bytes={"geno": 100, "snp": 200, "ind": 50, "pseudohaploid_json": 500},
        )
    text = buf.getvalue()
    assert "Stage 1" not in text
    assert "Stage 2" not in text
    assert "Stage 4" not in text
    assert "Stage 3" in text
    assert "Coverage report" in text
