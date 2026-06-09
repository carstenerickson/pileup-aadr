"""Shared dataclasses for the orchestrator.

`ExtractCliArgs` is a frozen view of the `extract` subcommand's CLI args; the
click decorator collects everything into a kwargs dict that maps directly into
this constructor. `ExtractCliArgs(**click_kwargs)` is the call site.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, get_args

# pileupCaller genotype-calling mode. Default randomHaploid: it pseudo-haploidizes
# the modern WGS target (one random read per site → 0% het) to MATCH the AADR
# 1240K panel, which is itself pseudo-haploid. majorityCall is also pseudo-haploid
# (consensus allele → 0% het). randomDiploid samples two reads → a DIPLOID, het-
# bearing call (~13% het on modern WGS) that does NOT match the panel — it is the
# legacy escape hatch only, and the .pseudohaploid sidecar records pseudohaploid=0
# when it is selected. See pileup_call.run_pileup_call.
#
# The Literal value of each mode is byte-identical to its pileupCaller CLI flag
# (`--{mode}`), which is what pileup_call relies on to build the argv — do not
# rename the members to e.g. snake_case without also remapping to the flag.
CallingMode = Literal["randomHaploid", "randomDiploid", "majorityCall"]

# Single source of truth for the selectable modes, in CLI-choice order. The click
# option and the `validate` flag-probe both derive from this, so adding a mode is
# a one-line edit to the Literal above.
CALLING_MODES: tuple[CallingMode, ...] = get_args(CallingMode)

# Modes whose output is pseudo-haploid (0% het, matches the AADR panel). This is an
# ALLOWLIST: any mode not listed here is treated as diploid (sidecar pseudohaploid=0,
# CLI warning fires), so a future calling mode fails CLOSED — it can never be silently
# mislabeled as pseudo-haploid, which was the original randomDiploid bug.
PSEUDOHAPLOID_MODES: frozenset[CallingMode] = frozenset({"randomHaploid", "majorityCall"})


def mode_is_pseudohaploid(mode: CallingMode) -> bool:
    """True if `mode`'s genotypes are pseudo-haploid (0% het) and match the AADR panel."""
    return mode in PSEUDOHAPLOID_MODES


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
    picard_shards: int | None = None  # None → memory-aware default at runtime
    min_mapq: int = 30
    min_baseq: int = 30
    no_baq: bool = False  # CLI flips --enable-baq → no_baq=False at the click layer
    seed: int = 42
    calling_mode: CallingMode = "randomHaploid"

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


__all__ = [
    "CALLING_MODES",
    "PSEUDOHAPLOID_MODES",
    "CallingMode",
    "CoverageCliArgs",
    "ExtractCliArgs",
    "mode_is_pseudohaploid",
]
