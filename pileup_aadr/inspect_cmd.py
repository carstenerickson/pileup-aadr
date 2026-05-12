"""`inspect` subcommand — click decorator. Implementation lives in `inspect_impl.py`."""
from __future__ import annotations

from pathlib import Path

import click


@click.command()
@click.argument(
    "aadr_snp", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "--json",
    "json_output",
    is_flag=True,
    help="Machine-readable JSON output (default: human-readable TSV)",
)
@click.pass_context
def inspect(ctx: click.Context, aadr_snp: Path, json_output: bool) -> None:
    """Structured summary of an AADR .snp file. Pure-Python; no external dependencies.

    AADR_SNP   AADR .snp file to inspect

    Fields emitted:
      total_rows                 int
      duplicate_rsid_count       int (always 0 — parser raises on duplicates)
      build                      "hg19" | "hg38" | "unknown"
      chrom_distribution         dict[chrom_str, count]
      allele_distribution        dict[ref_alt_pair, count]
      palindrome_count           int (A/T + C/G pairs)
      palindrome_fraction        float
      non_snp_count              int (always 0 — parser raises on non-ACGT)
      morgans_present            bool
      panel_guess                "1240k" | "HO" | "unknown"
    """
    from .inspect_impl import run_inspect

    exit_code = run_inspect(aadr_snp=aadr_snp, json_output=json_output)
    ctx.exit(exit_code)
