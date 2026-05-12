"""`validate` subcommand — click decorator. Implementation lives in `validate_impl.py`."""
from __future__ import annotations

from pathlib import Path

import click


@click.command()
@click.argument("bam", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument(
    "aadr_snp", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "--chain",
    "chain_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Chain file path (default: package-bundled)",
)
@click.option(
    "--ref-fasta",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Target FASTA matching BAM build",
)
@click.option(
    "--bam-build", type=click.Choice(["auto", "hg19", "hg38"]), default="auto"
)
@click.option(
    "--aadr-build", type=click.Choice(["auto", "hg19", "hg38"]), default="auto"
)
@click.option(
    "--json",
    "json_output",
    is_flag=True,
    help="Machine-readable JSON output",
)
@click.option(
    "-o",
    "--out",
    "output_prefix",
    type=click.Path(path_type=Path),
    default=None,
    help="Check this prefix for lock-held / overwrite-conflict (skip check if not given)",
)
@click.pass_context
def validate(
    ctx: click.Context,
    bam: Path,
    aadr_snp: Path,
    chain_path: Path | None,
    ref_fasta: Path | None,
    bam_build: str,
    aadr_build: str,
    json_output: bool,
    output_prefix: Path | None,
) -> None:
    """Pre-flight check (10 PASS/WARN/FAIL/SKIP checks).

    Does NOT run mpileup or pileupCaller. Same --chain / --ref-fasta / --bam-build /
    --aadr-build flags as `extract`.

    BAM       Aligned BAM/CRAM (hg19 or hg38, auto-detected from @SQ)
    AADR_SNP  AADR .snp file

    Output: TSV by default; --json for machine-readable.
    Exit 0 if all PASS; exit 1 if any FAIL.
    """
    from .validate_impl import run_validate

    exit_code = run_validate(
        bam=bam,
        aadr_snp=aadr_snp,
        chain_path=chain_path,
        ref_fasta=ref_fasta,
        bam_build_override=bam_build,  # type: ignore[arg-type]
        aadr_build_override=aadr_build,  # type: ignore[arg-type]
        json_output=json_output,
        output_prefix=output_prefix,
    )
    ctx.exit(exit_code)
