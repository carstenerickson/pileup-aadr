"""`extract` subcommand — thin click wrapper that delegates to `run_extract`.

The click decorator collects ~30 options into a kwargs dict that maps directly
into `ExtractCliArgs`. The orchestrator (`extract_orch.run_extract`) does the
real work.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import click

from .extract_orch import run_extract
from .types import ExtractCliArgs


@click.command()
@click.argument("bam", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument(
    "aadr_snp", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--out",
    "output_prefix",
    required=True,
    type=click.Path(path_type=Path),
    help="Output EIGENSTRAT prefix (.geno/.snp/.ind/.pseudohaploid.json written)",
)
@click.option(
    "--sample-name",
    type=str,
    default=None,
    help="Sample IID in output .ind (default: BAM @RG SM: -> filename stem)",
)
@click.option(
    "--pop",
    "pop_name",
    type=str,
    default=None,
    help="POP column in output .ind (default: --sample-name)",
)
@click.option(
    "--sex",
    type=click.Choice(["M", "F", "U"]),
    default="U",
    help="SEX column in output .ind",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing output files (default: refuse with exit 4)",
)
@click.option(
    "--chain",
    "chain_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Chain file path (default: package-bundled hg19ToHg38.over.chain.gz)",
)
@click.option(
    "--ref-fasta",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Target FASTA matching BAM build (default: auto-detect from BAM @PG)",
)
@click.option(
    "--bam-build",
    type=click.Choice(["auto", "hg19", "hg38"]),
    default="auto",
    help="Override BAM build detection",
)
@click.option(
    "--aadr-build",
    type=click.Choice(["auto", "hg19", "hg38"]),
    default="auto",
    help="Override AADR .snp build detection (default: auto; v0.1 AADR is hg19-native)",
)
@click.option(
    "--picard-mem",
    type=str,
    default="3g",
    help="Picard JVM heap size (default: 3g; verified empirically vs 2g floor)",
)
@click.option(
    "--picard-max-records",
    type=int,
    default=100_000,
    help="Picard MAX_RECORDS_IN_RAM (default: 100000; lower than Picard's 500000 default)",
)
@click.option(
    "--strict-chain-sha",
    is_flag=True,
    help="Enforce package-pinned SHA256 on user-supplied --chain (default: skip check)",
)
@click.option(
    "--insecure-chain",
    is_flag=True,
    help="Skip SHA256 verification entirely with stderr warning",
)
@click.option(
    "--keep-palindromes",
    is_flag=True,
    help="Pass A/T and C/G SNPs through Stage 1 (default: drop ~8%)",
)
@click.option(
    "--keep-alt-contigs",
    is_flag=True,
    help=(
        "Pass alt-haplotype contigs through Stage 2 "
        "(WARNING: pileupCaller will crash mid-Stage-3)"
    ),
)
@click.option(
    "--threads",
    type=click.IntRange(min=1),
    default=1,
    help=(
        "Number of parallel Stage 3 shards (default: 1 = single-process, "
        "byte-identical to v0.2; N>1 splits by chromosome and runs N "
        "concurrent samtools mpileup | pileupCaller processes)"
    ),
)
@click.option(
    "--picard-shards",
    type=click.IntRange(min=1),
    default=None,
    help=(
        "Stage 1 Picard concurrency. Default: derived from --threads and "
        "available memory (each Picard JVM is ~3 GB). Override only if you "
        "need to constrain or expand Stage 1 separately from Stage 3."
    ),
)
@click.option("--min-mapq", type=int, default=30, help="mpileup -q (default: 30)")
@click.option("--min-baseq", type=int, default=30, help="mpileup -Q (default: 30)")
@click.option(
    "--enable-baq",
    is_flag=True,
    help=(
        "Enable samtools BAQ by omitting the default -B flag (default: -B is passed, "
        "disabling samtools BAQ to match pileupCaller's recommended cmdline)"
    ),
)
@click.option(
    "--seed",
    type=int,
    default=42,
    help="pileupCaller --randomDiploid seed (default: 42)",
)
@click.option(
    "--liftover-yield-fail-pct",
    type=float,
    default=70.0,
    help="Exit 1 if Picard yield < N%% (default: 70)",
)
@click.option(
    "--liftover-yield-warn-pct",
    type=float,
    default=95.0,
    help="Stderr warning if yield < N%% (default: 95)",
)
@click.option(
    "--min-coverage",
    type=int,
    default=500_000,
    help="Exit 1 if non-missing autosomal calls < N (default: 500000)",
)
@click.option(
    "--warn-coverage",
    type=int,
    default=800_000,
    help="Stderr warning + JSON flag if calls < N (default: 800000)",
)
@click.option(
    "--report-json",
    type=click.Path(path_type=Path),
    default=None,
    help="Run-level summary JSON (consumed by ancestry-pipeline-tool gate)",
)
@click.option(
    "--report-tsv",
    type=click.Path(path_type=Path),
    default=None,
    help="Per-variant action TSV (streamed; constant memory)",
)
@click.option(
    "--tempdir",
    type=click.Path(path_type=Path),
    default=None,
    help="Override $TMPDIR for intermediates",
)
@click.option(
    "--keep-tempdir",
    is_flag=True,
    help="Always retain intermediates regardless of exit reason (default: retain only on crash)",
)
@click.option(
    "--clean-tempdir-on-crash",
    is_flag=True,
    help="Clean tempdir even on crash (for ephemeral CI/container environments)",
)
@click.pass_context
def extract(ctx: click.Context, **kwargs: Any) -> None:
    """Extract pseudohaploid genotypes at AADR sites from a BAM/CRAM.

    BAM       Aligned BAM/CRAM (hg19 or hg38, auto-detected from @SQ)
    AADR_SNP  AADR .snp file (hg19 coordinates through v66)

    The 4-stage pipeline:
      1. Lift AADR sites hg19 -> hg38 via Picard LiftoverVcf RECOVER_SWAPPED_REF_ALT
      2. Transform lifted VCF -> pileupCaller .snp + BED (with alt-contig filter)
      3. samtools mpileup | pileupCaller --randomDiploid
      4. Rejoin hg19 coordinates by rsID + invert dosage at SwappedAlleles

    For hg19-native BAMs, Stages 1/2/4 are skipped (no-lift fast path).
    """
    # CLI -> orchestrator BAQ flag mapping. The variable name `no_baq` is
    # awkwardly inverted (from the orchestrator's perspective): the value
    # `no_baq=True` means "DON'T pass `-B` to mpileup" — which actually
    # ENABLES samtools BAQ recalibration (-B is the disable-BAQ flag). So:
    #   --enable-baq passed (samtools BAQ wanted)  -> no_baq=True
    #   --enable-baq omitted (default; -B passed)  -> no_baq=False
    # --enable-baq passed → BAQ on → no_baq=False is wrong, flip: no_baq=True.
    # CLI default: BAQ off (samtools -B), so enable_baq=False → no_baq=False.
    enable_baq = kwargs.pop("enable_baq", False)
    kwargs["no_baq"] = enable_baq
    args = ExtractCliArgs(**kwargs)
    exit_code = run_extract(args)
    ctx.exit(exit_code)
