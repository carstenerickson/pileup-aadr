"""Shared dataclasses for the orchestrator.

`ExtractCliArgs` is a frozen view of the `extract` subcommand's CLI args; the
click decorator collects everything into a kwargs dict that maps directly into
this constructor. `ExtractCliArgs(**click_kwargs)` is the call site.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal


@dataclass(frozen=True)
class ExtractCliArgs:
    """Frozen view of the `extract` subcommand's CLI args; mirrors click 1:1."""

    # Required positionals
    bam: Path
    aadr_snp: Path

    # Output
    output_prefix: Path
    sample_name: str | None = None
    pop_name: str | None = None
    sex: Literal["M", "F", "U"] = "U"
    overwrite: bool = False

    # Liftover
    chain_path: Path | None = None
    ref_fasta: Path | None = None
    bam_build: Literal["auto", "hg19", "hg38"] = "auto"
    aadr_build: Literal["auto", "hg19", "hg38"] = "auto"
    picard_mem: str = "3g"
    picard_max_records: int = 100_000
    strict_chain_sha: bool = False
    insecure_chain: bool = False

    # Filtering
    keep_palindromes: bool = False
    keep_alt_contigs: bool = False

    # Pileup / call
    threads: int = 1
    min_mapq: int = 30
    min_baseq: int = 30
    no_baq: bool = False  # CLI flips --enable-baq → no_baq=False at the click layer
    seed: int = 42

    # Validation thresholds
    liftover_yield_fail_pct: float = 70.0
    liftover_yield_warn_pct: float = 95.0
    min_coverage: int = 500_000
    warn_coverage: int = 800_000

    # Reporting
    report_json: Path | None = None
    report_tsv: Path | None = None

    # Tempdir
    tempdir: Path | None = None
    keep_tempdir: bool = False
    clean_tempdir_on_crash: bool = False


@dataclass(frozen=True)
class CoverageCliArgs:
    """Frozen view of the `coverage` subcommand's CLI args."""

    bam: Path
    regions: Path | None = None
    threads: int = 4
    quantize: str = "0:1:5:10:30"
    json_output: bool = False


__all__ = ["CoverageCliArgs", "ExtractCliArgs"]
