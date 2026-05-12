"""CLI entry point.

Per HLD §"Subcommand structure" + LLD §3:
    pileup-aadr extract   BAM AADR_SNP -o PREFIX [...]   the canonical 4-stage pipeline
    pileup-aadr validate  BAM AADR_SNP                   pre-flight check
    pileup-aadr coverage  BAM [--regions BED]            per-chrom coverage report
    pileup-aadr inspect   AADR_SNP                       structured summary of an AADR .snp
"""
from __future__ import annotations

import logging
import sys
import traceback

import click

from . import __version__
from .errors import PileupAadrError, format_error
from .extract_cmd import extract
from .inspect_cmd import inspect
from .logging_config import configure_logging
from .validate_cmd import validate

# `coverage` subcommand deferred to Days 6-7 per the project plan; placeholder click stub
# kept off the root group until then.


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="pileup-aadr")
@click.option("--verbose", is_flag=True, help="DEBUG-level logging to stderr")
@click.option(
    "--quiet",
    is_flag=True,
    help="Suppress INFO messages; warnings still go to stderr",
)
@click.pass_context
def cli(ctx: click.Context, verbose: bool, quiet: bool) -> None:
    """Extract pseudohaploid genotypes at AADR sites from a user BAM."""
    if verbose and quiet:
        raise click.UsageError("--verbose and --quiet are mutually exclusive")
    level = logging.DEBUG if verbose else (logging.WARNING if quiet else logging.INFO)
    configure_logging(level=level)
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet


cli.add_command(extract)
cli.add_command(validate)
cli.add_command(inspect)


def main() -> int:
    """Entry point (per pyproject.toml [project.scripts]). Wraps cli() with error formatter.

    Exit codes (per HLD §"Exit codes"):
        0  success
        1  soft-validation failure (LiftoverYieldError, CoverageGateFailure)
        2  I/O failure (chain/FASTA/BAM not found, subprocess crashed, lock held)
        3  invariant violation (build mismatch, AADR malformed, defensive sanity check)
        4  usage error (bad CLI args, missing/wrong-version external binary)
    """
    try:
        cli.main(standalone_mode=False)
        return 0
    except click.exceptions.Exit as e:
        return e.exit_code
    except click.UsageError as e:
        e.show()
        return e.exit_code or 4
    except click.ClickException as e:
        e.show()
        return e.exit_code or 4
    except PileupAadrError as e:
        sys.stderr.write(format_error(e))
        return e.exit_code
    except Exception as e:
        sys.stderr.write(f"pileup-aadr: [crash] {type(e).__name__}: {e}\n")
        traceback.print_exc(file=sys.stderr)
        return 3  # treat uncaught as invariant violation


if __name__ == "__main__":
    sys.exit(main())
