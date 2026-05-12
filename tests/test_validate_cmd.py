"""Smoke tests for the `validate` subcommand.

Focus on the file-side checks (BAM parse, AADR parse, output prefix) since the
binary-version probes are stubbed for Day 1 (full impl lands on Day 2 with tool_wrapper.py).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pileup_aadr.cli import cli


# Minimal fake BAM: just gzip+BAM-magic. Real validation requires an actual BAM with
# @SQ chr1; for these tests we write one programmatically.
@pytest.fixture
def minimal_bam(tmp_path: Path) -> Path:
    """Construct a real (but tiny) BAM file with @SQ chr1 + a single read.

    Using pysam directly so the file passes BAM-format detection and has a parseable
    header. We don't need any reads — header alone exercises the validate pre-flight.
    """
    import pysam

    out = tmp_path / "tiny.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        # chr1 length matches hg38 within ±1Mb so detect_bam_build returns "hg38"
        "SQ": [{"SN": "chr1", "LN": 248_956_422}],
        "RG": [{"ID": "rg1", "SM": "test_sample", "LB": "lib", "PL": "ILLUMINA"}],
    }
    with pysam.AlignmentFile(str(out), "wb", header=header):
        pass  # zero reads
    pysam.sort("-o", str(out), str(out))
    pysam.index(str(out))
    return out


def test_validate_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["validate", "--help"])
    assert result.exit_code == 0
    assert "BAM" in result.output
    assert "AADR_SNP" in result.output


def test_validate_tsv_runs(minimal_bam: Path, aadr_chr22_slice: Path) -> None:
    """Default TSV output runs end-to-end and produces a result table."""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["--quiet", "validate", str(minimal_bam), str(aadr_chr22_slice), "--bam-build", "hg38"],
    )
    # Exit code may be 0 (all PASS), 1 (some FAIL — likely because samtools/picard/etc.
    # may or may not be on PATH in the test env). The important thing is the command ran
    # without crashing and produced a parseable TSV.
    assert result.exit_code in (0, 1), result.output
    lines = result.output.strip().splitlines()
    assert lines[0] == "status\tcheck\tdetail"
    statuses = {line.split("\t", 1)[0] for line in lines[1:]}
    # Every entry must be one of PASS/WARN/FAIL/SKIP
    assert statuses <= {"PASS", "WARN", "FAIL", "SKIP"}


def test_validate_json_output(minimal_bam: Path, aadr_chr22_slice: Path) -> None:
    """--json emits a parseable object with checks + exit_code."""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--quiet",
            "validate",
            "--json",
            str(minimal_bam),
            str(aadr_chr22_slice),
            "--bam-build",
            "hg38",
        ],
    )
    assert result.exit_code in (0, 1), result.output
    data = json.loads(result.output)
    assert "checks" in data
    assert "exit_code" in data
    assert isinstance(data["checks"], list)
    for check in data["checks"]:
        assert set(check) == {"name", "status", "detail"}
        assert check["status"] in {"PASS", "WARN", "FAIL", "SKIP"}


def test_validate_aadr_check_passes(
    minimal_bam: Path, aadr_chr22_slice: Path
) -> None:
    """The AADR slice parses cleanly so its parse check should be PASS."""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--quiet",
            "validate",
            "--json",
            str(minimal_bam),
            str(aadr_chr22_slice),
            "--bam-build",
            "hg38",
        ],
    )
    data = json.loads(result.output)
    aadr_check = next(c for c in data["checks"] if c["name"] == "AADR .snp parseable")
    assert aadr_check["status"] == "PASS"
    assert "50 unique rsIDs" in aadr_check["detail"]


def test_validate_no_lift_skips_picard_java(
    minimal_bam: Path, tmp_path: Path
) -> None:
    """When AADR_build == BAM_build (both hg19 here), picard + java checks SKIP."""
    # Make a tiny hg19-detectable AADR
    snp = tmp_path / "hg19_aadr.snp"
    snp.write_text("rs1 1 0.0 249200000 A G\n")  # closer to hg19 chr1 end than hg38
    # Construct a hg19-length BAM
    import pysam

    bam = tmp_path / "hg19.bam"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": 249_250_621}],  # hg19 length
    }
    with pysam.AlignmentFile(str(bam), "wb", header=header):
        pass
    pysam.sort("-o", str(bam), str(bam))
    pysam.index(str(bam))

    runner = CliRunner()
    result = runner.invoke(
        cli, ["--quiet", "validate", "--json", str(bam), str(snp), "--bam-build", "hg19", "--aadr-build", "hg19"]
    )
    data = json.loads(result.output)
    picard_check = next(c for c in data["checks"] if c["name"] == "picard binary")
    java_check = next(c for c in data["checks"] if c["name"] == "java binary")
    assert picard_check["status"] == "SKIP"
    assert java_check["status"] == "SKIP"
    assert "no-lift" in picard_check["detail"]
