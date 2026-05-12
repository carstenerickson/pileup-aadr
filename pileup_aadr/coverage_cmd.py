"""`coverage` subcommand — thin click wrapper around `coverage_impl.run_coverage`."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import click

from .coverage_impl import run_coverage
from .types import CoverageCliArgs


@click.command()
@click.argument("bam", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--regions",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="BED file restricting coverage to specific regions (e.g., AADR positions)",
)
@click.option(
    "--threads",
    type=int,
    default=4,
    help="mosdepth worker threads (default: 4)",
)
@click.option(
    "--quantize",
    type=str,
    default="0:1:5:10:30",
    help=(
        "mosdepth quantize bins (default: 0:1:5:10:30 — split coverage into "
        "0/1-4/5-9/10-29/30+ buckets)"
    ),
)
@click.option(
    "--json",
    "json_output",
    is_flag=True,
    help="Emit structured JSON instead of TSV",
)
@click.pass_context
def coverage(ctx: click.Context, **kwargs: Any) -> None:
    """Per-chromosome BAM coverage report via mosdepth.

    BAM   Aligned BAM/CRAM file (must be indexed)

    Diagnostic output reached for when `pileup-aadr extract`'s coverage gate
    FAILs — distinguishes "BAM is low-coverage globally" from "BAM is fine
    but the AADR 1240k positions specifically have poor coverage".
    """
    args = CoverageCliArgs(**kwargs)
    exit_code = run_coverage(args)
    ctx.exit(exit_code)
