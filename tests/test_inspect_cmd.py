"""Smoke tests for the `inspect` subcommand end-to-end via CliRunner."""
from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from pileup_aadr.cli import cli


def test_inspect_tsv_output(aadr_chr22_slice: Path) -> None:
    """Default TSV output emits all 10 fields."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--quiet", "inspect", str(aadr_chr22_slice)])
    assert result.exit_code == 0, result.output
    lines = result.output.strip().splitlines()
    assert lines[0] == "field\tvalue"
    fields = {line.split("\t", 1)[0] for line in lines[1:]}
    expected_fields = {
        "total_rows",
        "duplicate_rsid_count",
        "build",
        "chrom_distribution",
        "chrom_set",
        "allele_distribution",
        "palindrome_count",
        "palindrome_fraction",
        "non_snp_count",
        "morgans_present",
        "panel_guess",
    }
    assert fields == expected_fields


def test_inspect_json_output(aadr_chr22_slice: Path) -> None:
    """--json emits a parseable object with the expected shape."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--quiet", "inspect", "--json", str(aadr_chr22_slice)])
    assert result.exit_code == 0, result.output
    data = json.loads(result.output)
    assert data["total_rows"] == 50
    assert data["duplicate_rsid_count"] == 0
    # chr22 slice — should map to "22" in the chrom_distribution
    assert data["chrom_distribution"] == {"22": 50}
    # The slice has no chr1 rows so build detection raises and we fall back to "unknown"
    assert data["build"] == "unknown"
    # Slice is 50 rows — too small for either 1240k or HO heuristic
    assert data["panel_guess"] == "unknown"
    assert data["non_snp_count"] == 0
    # All gen_morgans values in the slice are 0.0, so morgans_present is False
    assert data["morgans_present"] is False
    # chr22 only → not autosomes_only (which requires all 22), not chrY/chrM only → custom
    assert data["chrom_set"] == "custom"


def test_inspect_palindrome_count(tmp_path: Path) -> None:
    """A panel with explicit A/T + C/G alleles flags the palindrome count correctly."""
    snp = tmp_path / "palindromes.snp"
    snp.write_text(
        "rs1 22 0.0 1000 A T\n"  # A/T palindrome
        "rs2 22 0.0 2000 T A\n"  # T/A palindrome
        "rs3 22 0.0 3000 C G\n"  # C/G palindrome
        "rs4 22 0.0 4000 G C\n"  # G/C palindrome
        "rs5 22 0.0 5000 A G\n"  # not a palindrome
        "rs6 22 0.0 6000 C T\n"  # not a palindrome
    )
    runner = CliRunner()
    result = runner.invoke(cli, ["--quiet", "inspect", "--json", str(snp)])
    assert result.exit_code == 0, result.output
    data = json.loads(result.output)
    assert data["palindrome_count"] == 4
    assert data["palindrome_fraction"] == round(4 / 6, 4)


def test_inspect_help() -> None:
    """`inspect --help` exits 0 and mentions the AADR_SNP positional."""
    runner = CliRunner()
    result = runner.invoke(cli, ["inspect", "--help"])
    assert result.exit_code == 0
    assert "AADR_SNP" in result.output
