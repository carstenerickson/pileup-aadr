"""Counter dataclasses populated by each stage and aggregated by the orchestrator.

Shape mirrors the v2.1 JSON report schema 1:1 — `ExtractCounters` serializes via
`dataclasses.asdict()` directly into the schema-versioned JSON output, no translation
layer.

Frozen dataclasses (the leaf counters) are constructed once and never mutated.
`ExtractCounters` itself is mutable so the orchestrator can replace its `coverage`
and `gates` fields after the gate-evaluation pass without rebuilding the whole tree.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Stage1InputFilters:
    """Counter return shape from `sites_vcf.build_sites_vcf`, fed into Stage1LiftCounters
    as the `input_filters` dict-equivalent fields."""

    palindrome_drops: int
    non_snp_drops: int
    non_autosome_drops: int  # 0 by default (we keep X/Y); reserved for autosome-only mode
    rows_written: int


@dataclass(frozen=True)
class Stage1LiftCounters:
    """Populated by `lift.lift_aadr_sites()`."""

    wallclock_seconds: float
    input_sites_after_filters: int  # denominator for liftover_yield_pct (post-input-filter)
    lifted_sites: int  # numerator (Picard's OUTPUT record count)
    liftover_yield_pct: float  # 100 * lifted_sites / input_sites_after_filters
    liftover_yield_warning: bool  # True if yield < --liftover-yield-warn-pct
    rejected_by_reason: dict[str, int]
    # ^ {"NoTarget": N, "MismatchedRefAllele": N, "IndelStraddlesMultipleIntervals": N,
    #    "SwappedAlleles": N, "other": N}
    swapped_alleles_count: int  # Picard's SwappedAlleles INFO marker count (recovered, not rejected)
    input_filters: dict[str, int]
    # ^ {"palindrome_drops": N, "non_snp_drops": N, "non_autosome_drops": N}


@dataclass(frozen=True)
class Stage2TransformCounters:
    """Populated by `transform.build_pileupcaller_snp_and_bed()`."""

    wallclock_seconds: float
    alt_contig_drops: int
    output_sites: int


@dataclass(frozen=True)
class PileupCallerSummary:
    """Parsed verbatim from pileupCaller's stderr summary stats TSV (v2.1 verified).

    Field names match pileupCaller's stderr column headers exactly so the parser
    can populate via splat.
    """

    total_sites: int
    non_missing_calls: int
    avg_raw_reads: float
    avg_damage_cleaned_reads: float
    avg_sampled_from: float


@dataclass(frozen=True)
class Stage3CallCounters:
    """Populated by `pileup_call.run_pileup_call()`."""

    wallclock_seconds: float
    pileupcaller_summary: PileupCallerSummary


@dataclass(frozen=True)
class Stage4RejoinCounters:
    """Populated by `rejoin.rejoin_aadr_frame()` (or `_no_lift_fast_path_finalize` for fast path).

    Invariant: `ref_alt_swap_count` == `Stage1LiftCounters.swapped_alleles_count` for matched
    rsIDs in the lift-path case (cross-checked at orchestrator level for diagnostic; mismatch
    is surfaced as a WARNING but doesn't fail the run).
    """

    wallclock_seconds: float
    rsid_matched: int
    ref_alt_swap_count: int
    allele_mismatch_drops: int  # defensive sanity check; ~0 expected
    output_variants: int


@dataclass(frozen=True)
class CoverageCounters:
    """Populated by rejoin (lift path) or `_no_lift_fast_path_finalize` (fast path)."""

    non_missing_autosomal_calls: int  # gated count: chr1-22 only
    coverage_fraction: float  # diagnostic: non_missing_autosomal_calls / autosomal_aadr_count
    coverage_warning: bool  # True if non_missing_autosomal_calls < --warn-coverage
    per_chrom_call_count: dict[str, int]  # all chroms present in output (chr1-22, chrX, chrY, chrM)


@dataclass
class ExtractCounters:
    """Aggregated by the orchestrator (`extract_orch.run_extract`).

    Not frozen — orchestrator replaces `coverage` and `gates` after the gate-evaluation
    pass.

    For the no-lift fast path (AADR build == BAM build), Stages 1/2/4 are skipped and
    their counters are None. The orchestrator's JSON serialization omits null stages
    cleanly via a custom default that drops None-valued top-level fields.
    """

    stage_1_lift: Stage1LiftCounters | None
    stage_2_transform: Stage2TransformCounters | None
    stage_3_call: Stage3CallCounters
    stage_4_rejoin: Stage4RejoinCounters | None
    coverage: CoverageCounters
    gates: dict[str, str]  # {"liftover_yield": "PASS"|"WARN"|"FAIL"|"N/A", "coverage": same}
    wallclock_total_seconds: float


__all__ = [
    "CoverageCounters",
    "ExtractCounters",
    "PileupCallerSummary",
    "Stage1InputFilters",
    "Stage1LiftCounters",
    "Stage2TransformCounters",
    "Stage3CallCounters",
    "Stage4RejoinCounters",
]
