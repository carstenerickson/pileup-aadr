"""Stage 1: Picard LiftoverVcf + chain-file resolution.

Two layers:

1. Chain-file resolution (Day 2): `get_bundled_chain_path` + `chain_file_path` +
   `resolve_chain_for_extract`. Bundled chain is SHA-verified at startup.

2. Stage 1 lift (Day 3): `lift_aadr_sites` runs Picard via `ToolWrapper`, parses
   Picard's structured stderr + the `_rejected.vcf` file to populate
   `Stage1LiftCounters`, evaluates the liftover-yield gate, and produces the
   lifted VCF that Stages 2-4 consume.
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import importlib.resources as res
import logging
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

import pysam

from .counters import Stage1InputFilters, Stage1LiftCounters, Stage2TransformCounters
from .errors import (
    ChainFileNotFound,
    ChainFileSHAError,
    LiftoverYieldError,
    PileupAadrInternalError,
)
from .tool_wrapper import PICARD_SPEC, ToolWrapper

log = logging.getLogger(__name__)


def get_bundled_chain_path() -> Path:
    """Resolve the bundled chain file path; verify SHA at access.

    Always-on verification catches "wheel got corrupted in transit" at ~2 ms cost
    (read 223 KB + SHA-256). Per HLD §"Bundled chain file packaging":

      - sha256 sidecar lives at pileup_aadr/data/hg19ToHg38.over.chain.gz.sha256
      - Mismatch → ChainFileSHAError + reinstall guidance

    Returns:
        Path to the bundled chain file (materialized via importlib.resources.as_file
        if the package is in a zip archive; otherwise direct filesystem path).

    Raises:
        ChainFileSHAError: bundled chain bytes don't match the .sha256 sidecar
            (corrupt install).
    """
    chain_resource = res.files("pileup_aadr") / "data" / "hg19ToHg38.over.chain.gz"
    sha_resource = res.files("pileup_aadr") / "data" / "hg19ToHg38.over.chain.gz.sha256"

    expected_sha = sha_resource.read_text().strip().split()[0]
    chain_bytes = chain_resource.read_bytes()
    actual_sha = hashlib.sha256(chain_bytes).hexdigest()

    if actual_sha != expected_sha:
        raise ChainFileSHAError(
            what=f"bundled chain {chain_resource}",
            why=f"SHA-256 mismatch: expected {expected_sha[:16]}…, got {actual_sha[:16]}…",
            fix="Reinstall pileup-aadr (the wheel may have been corrupted in transit)",
        )

    # For filesystem-installed packages, files() returns a regular Path. For zip-archived
    # packages we'd need res.as_file() to materialize a real path; for v0.1 we assume
    # filesystem install (the common case for pip install -e and wheel installs).
    return Path(str(chain_resource))


def chain_file_path(
    cli_chain: Path | None,
    env_chain_dir: Path | None,
    *,
    strict_sha: bool = False,
    insecure: bool = False,
) -> Path:
    """3-tier chain-file resolution per HLD §"Chain & reference dependencies".

    Resolution order:
      1. `cli_chain` (--chain PATH explicit override)
      2. `env_chain_dir` / hg19ToHg38.over.chain.gz ($PILEUP_AADR_CHAIN_DIR env var)
      3. Bundled chain at pileup_aadr/data/ (always present after install; SHA-verified)

    Args:
        cli_chain: --chain PATH if user passed it
        env_chain_dir: $PILEUP_AADR_CHAIN_DIR resolved Path (None if not set)
        strict_sha: if True, run the bundled-chain SHA check on user-supplied chains too.
            Default False — explicit user choice is trusted unless --strict-chain-sha set.
        insecure: if True, skip SHA verification entirely with a stderr WARNING.
            Default False — bundled-chain SHA is always verified for safety.

    Returns:
        Path to the chain file to use. The bundled path returned has already been
        SHA-verified; user-supplied paths are returned without verification unless
        strict_sha is True.

    Raises:
        ChainFileNotFound: cli_chain was given but the path doesn't exist.
        ChainFileSHAError: bundled-chain SHA mismatch, or strict_sha + user-chain mismatch.
    """
    if cli_chain is not None:
        if not cli_chain.exists():
            raise ChainFileNotFound(
                what=str(cli_chain),
                why="--chain path does not exist",
                fix="Check the path; or omit --chain to use the bundled hg19ToHg38.over.chain.gz",
            )
        if strict_sha and not insecure:
            _verify_user_chain_sha(cli_chain)
        elif insecure:
            log.warning(
                "Chain SHA verification skipped per --insecure-chain (using %s)",
                cli_chain,
            )
        return cli_chain

    if env_chain_dir is not None:
        candidate = env_chain_dir / "hg19ToHg38.over.chain.gz"
        if candidate.exists():
            log.debug("Using chain from $PILEUP_AADR_CHAIN_DIR: %s", candidate)
            return candidate

    return get_bundled_chain_path()


def _verify_user_chain_sha(user_chain: Path) -> None:
    """Verify a user-supplied --chain PATH matches the package-pinned SHA.

    Used only when --strict-chain-sha is set. The pinned SHA is read from the
    bundled .sha256 sidecar; the user-supplied chain is hashed and compared.
    """
    sha_resource = res.files("pileup_aadr") / "data" / "hg19ToHg38.over.chain.gz.sha256"
    expected_sha = sha_resource.read_text().strip().split()[0]
    actual_sha = hashlib.sha256(user_chain.read_bytes()).hexdigest()
    if actual_sha != expected_sha:
        raise ChainFileSHAError(
            what=str(user_chain),
            why=(
                f"--strict-chain-sha enforced and user-supplied chain SHA differs from "
                f"the package-pinned canonical SHA: expected {expected_sha[:16]}…, "
                f"got {actual_sha[:16]}…"
            ),
            fix=(
                "Either supply the canonical UCSC chain (matches the bundled SHA) or "
                "drop --strict-chain-sha to trust your --chain choice; or pass "
                "--insecure-chain to skip verification with a stderr warning"
            ),
        )


def resolve_chain_for_extract(
    cli_chain: Path | None,
    *,
    strict_sha: bool = False,
    insecure: bool = False,
) -> Path:
    """Convenience wrapper that pulls $PILEUP_AADR_CHAIN_DIR from env + delegates.

    Used by extract_orch (Day 5) and validate (Day 2) so the env-var lookup logic
    lives in one place.
    """
    env_dir = os.environ.get("PILEUP_AADR_CHAIN_DIR")
    return chain_file_path(
        cli_chain=cli_chain,
        env_chain_dir=Path(env_dir) if env_dir else None,
        strict_sha=strict_sha,
        insecure=insecure,
    )


# ---------------------------------------------------------------------------
# Stage 1: Picard LiftoverVcf
# ---------------------------------------------------------------------------


# Patterns verified empirically against Picard 3.3.0 stderr captures.
# v2.2 critique M14: split into REQUIRED (raise on missing) vs OPTIONAL (default
# on missing). The "swapped" line could plausibly be omitted by a future Picard
# version on a 0-swap run; defaulting to "0" keeps the parser robust without
# sacrificing fidelity (the rejected.vcf parser cross-validates).
_PICARD_REQUIRED_PATTERNS: Final[dict[str, re.Pattern[str]]] = {
    "processed": re.compile(r"Processed (\d+) variants\."),
    "failed": re.compile(r"(\d+) variants failed to liftover\."),
    "mismatched_ref": re.compile(
        r"(\d+) variants lifted over but had mismatching reference alleles after lift over\."
    ),
    "yield_pct_picard": re.compile(
        r"([\d.]+)% of variants were not successfully lifted over"
    ),
}
_PICARD_OPTIONAL_PATTERNS: Final[dict[str, tuple[re.Pattern[str], str]]] = {
    "swapped": (
        re.compile(r"(\d+) variants were lifted by swapping REF/ALT alleles\."),
        "0",
    ),
}

# Picard rejection FILTER values we recognize. SwappedAlleles is technically
# emitted on RECOVERED sites (in OUTPUT, not REJECT) but defensively bucketed
# here in case a future Picard version writes it to REJECT for some path.
_REJECT_FILTER_BUCKETS: Final[set[str]] = {
    "NoTarget",
    "MismatchedRefAllele",
    "IndelStraddlesMultipleIntervals",
    "SwappedAlleles",
}


def lift_aadr_sites(
    sites_vcf_path: Path,
    chain_path: Path,
    target_fasta_path: Path,
    output_lifted_vcf: Path,
    output_rejected_vcf: Path,
    input_filter_counters: Stage1InputFilters,
    *,
    picard_mem: str = "3g",
    picard_max_records: int = 100_000,
    yield_fail_pct: float = 70.0,
    yield_warn_pct: float = 95.0,
) -> Stage1LiftCounters:
    """Run Picard LiftoverVcf with REF_AWARE handling.

    Args:
        sites_vcf_path: input VCF (built by `sites_vcf.build_sites_vcf`).
        chain_path: resolved chain file (bundled or user-supplied).
        target_fasta_path: TARGET assembly FASTA (typically hg38). The .dict
            sidecar must exist alongside (Picard requires it).
        output_lifted_vcf: where to write Picard's OUTPUT.
        output_rejected_vcf: where to write Picard's REJECT.
        input_filter_counters: from `sites_vcf.build_sites_vcf` — fed into the
            returned counter (the input_filters dict) and used as the denominator
            for liftover_yield_pct.
        picard_mem: JVM heap (passed as -Xmx<value>). Default 3g per HLD §"Memory model".
        picard_max_records: --MAX_RECORDS_IN_RAM. Default 100K (Picard's 500K
            spills excessively for a 1.2M-site VCF).
        yield_fail_pct: exit 1 below this (default 70.0).
        yield_warn_pct: stderr WARNING + JSON flag below this (default 95.0).

    Returns:
        Stage1LiftCounters with parsed counts + wallclock + yield-gate result.

    Raises:
        LiftoverYieldError: yield < `yield_fail_pct`.
        ToolSubprocessError: Picard exited non-zero.
        PileupAadrInternalError: Picard's stderr couldn't be parsed
            (suggests Picard format change; pin tool version).
    """
    output_lifted_vcf.parent.mkdir(parents=True, exist_ok=True)
    output_rejected_vcf.parent.mkdir(parents=True, exist_ok=True)

    # H9 fix: jvm_args go via the dedicated parameter, NOT in args. ToolWrapper
    # injects them between `java` and `-jar` so we don't double-jar.
    jvm_args = [f"-Xmx{picard_mem}"]
    args = [
        "LiftoverVcf",
        "--INPUT", str(sites_vcf_path),
        "--OUTPUT", str(output_lifted_vcf),
        "--CHAIN", str(chain_path),
        "--REFERENCE_SEQUENCE", str(target_fasta_path),
        "--REJECT", str(output_rejected_vcf),
        "--RECOVER_SWAPPED_REF_ALT", "true",
        "--WRITE_ORIGINAL_POSITION", "true",
        "--WARN_ON_MISSING_CONTIG", "true",
        "--MAX_RECORDS_IN_RAM", str(picard_max_records),
    ]
    stderr_path = output_lifted_vcf.parent / "picard.stderr"

    wrapper = ToolWrapper(PICARD_SPEC)
    t0 = time.perf_counter()
    wrapper.run(
        args=args,
        jvm_args=jvm_args,
        capture_stderr_to=stderr_path,
        check=True,  # raises ToolSubprocessError on non-zero exit
    )
    wallclock = time.perf_counter() - t0

    # Parse Picard's structured stderr
    stderr_text = stderr_path.read_text()
    parsed = parse_picard_stderr(stderr_text)

    # Parse _rejected.vcf for FILTER categorization
    rejected_by_reason = parse_rejected_vcf(output_rejected_vcf)

    # Compute yield (denominator = post-input-filter rows the orchestrator wrote
    # to sites_vcf, NOT raw AADR row count).
    denom = input_filter_counters.rows_written
    numer = parsed["lifted_count"]
    yield_pct = 100.0 * numer / denom if denom > 0 else 0.0
    yield_warning = yield_pct < yield_warn_pct

    if yield_pct < yield_fail_pct:
        # Find dominant rejection reason for the diagnostic
        dominant = max(
            rejected_by_reason.items(),
            key=lambda kv: kv[1],
            default=("unknown", 0),
        )
        raise LiftoverYieldError(
            what=f"yield {yield_pct:.2f}% (< {yield_fail_pct}% gate)",
            why=(
                f"{numer}/{denom} sites lifted; dominant rejection: "
                f"{dominant[0]} ({dominant[1]} sites)"
            ),
            fix=(
                "Verify chain file matches BAM build + verify --ref-fasta matches "
                f"BAM @PG. For dominant rejection '{dominant[0]}', see the "
                "troubleshooting guide."
            ),
        )

    if yield_warning:
        dominant_name = max(
            rejected_by_reason.items(), key=lambda kv: kv[1], default=("none", 0)
        )[0]
        log.warning(
            "Liftover yield %.2f%% is below warn threshold %.1f%% (still proceeding); "
            "dominant rejection reason: %s",
            yield_pct,
            yield_warn_pct,
            dominant_name,
        )

    log.info(
        "Stage 1 complete: %d lifted of %d input (%.2f%% yield, %d swapped); wallclock %.1fs",
        numer,
        denom,
        yield_pct,
        parsed["swapped_count"],
        wallclock,
    )

    return Stage1LiftCounters(
        wallclock_seconds=wallclock,
        input_sites_after_filters=denom,
        lifted_sites=numer,
        liftover_yield_pct=round(yield_pct, 4),
        liftover_yield_warning=yield_warning,
        rejected_by_reason=rejected_by_reason,
        swapped_alleles_count=parsed["swapped_count"],
        input_filters=input_filter_counters,
    )


@dataclass(frozen=True)
class PicardShardSpec:
    """One Stage 1 shard — N shards per run, grouped by source chrom for balance."""

    shard_index: int
    source_chroms: tuple[str, ...]  # canonical "chr1".."chr22" / "chrX" / "chrY"
    input_vcf: Path               # per-shard sites VCF (full header + filtered body)
    lifted_vcf: Path              # per-shard Picard OUTPUT
    rejected_vcf: Path            # per-shard Picard REJECT
    stderr_path: Path             # per-shard Picard stderr


def build_picard_shard_manifest(
    sites_vcf_path: Path,
    shard_dir: Path,
    n_shards: int,
) -> list[PicardShardSpec]:
    """Partition sites VCF into N shards via LPT bin-packing on per-chrom counts.

    Two-pass: pass 1 counts records per source chrom; pass 2 streams records to
    per-shard input VCF files. All shards get the full VCF header (Picard validates
    each CHROM against declared ##contig lines). Clamps n_shards to chrom count.

    Args:
        sites_vcf_path: sorted sites VCF from build_sites_vcf.
        shard_dir: parent directory for shard_{00..NN}/ subdirs.
        n_shards: target shard count; clamped to number of chroms with sites.

    Returns:
        list[PicardShardSpec] sorted by shard_index, length = min(n_shards, chroms_with_sites).
    """
    # Pass 1: collect header lines + count records per source chrom
    header_lines: list[str] = []
    chrom_counts: dict[str, int] = {}
    with open(sites_vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                chrom = line.split("\t", 1)[0]
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

    if not chrom_counts:
        raise PileupAadrInternalError(
            what="build_picard_shard_manifest",
            why="sites VCF has no data records (build_sites_vcf produced empty output)",
            fix="Inspect Stage 1 input filters; AADR panel may have no canonical-chrom sites",
        )

    # Clamp to chrom count
    actual_shards = min(n_shards, len(chrom_counts))
    if actual_shards < n_shards:
        log.info(
            "Picard shard count clamped from %d to %d (only %d chroms with sites)",
            n_shards, actual_shards, len(chrom_counts),
        )

    # LPT bin-packing: sort chroms by count descending; greedily assign to least-loaded shard
    sorted_chroms = sorted(chrom_counts.keys(), key=lambda c: chrom_counts[c], reverse=True)
    shard_totals = [0] * actual_shards
    shard_chroms: list[list[str]] = [[] for _ in range(actual_shards)]
    for chrom in sorted_chroms:
        min_shard = min(range(actual_shards), key=lambda i: shard_totals[i])
        shard_chroms[min_shard].append(chrom)
        shard_totals[min_shard] += chrom_counts[chrom]

    chrom_to_shard: dict[str, int] = {
        c: i for i, chroms in enumerate(shard_chroms) for c in chroms
    }

    # Open per-shard input VCF files and write full header to each
    shard_dir.mkdir(parents=True, exist_ok=True)
    specs: list[PicardShardSpec] = []
    shard_handles: dict[int, Any] = {}
    try:
        for i, chroms in enumerate(shard_chroms):
            if not chroms:
                continue
            sd = shard_dir / f"shard_{i:02d}"
            sd.mkdir(exist_ok=True)
            input_vcf = sd / "input.vcf"
            specs.append(PicardShardSpec(
                shard_index=i,
                source_chroms=tuple(sorted(chroms)),
                input_vcf=input_vcf,
                lifted_vcf=sd / "lifted.vcf",
                rejected_vcf=sd / "rejected.vcf",
                stderr_path=sd / "picard.stderr",
            ))
            fh = open(input_vcf, "w")
            shard_handles[i] = fh
            for line in header_lines:
                fh.write(line)

        # Pass 2: stream body records to per-shard files
        with open(sites_vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                chrom = line.split("\t", 1)[0]
                shard_idx = chrom_to_shard.get(chrom)
                if shard_idx is not None:
                    shard_handles[shard_idx].write(line)
    finally:
        for fh in shard_handles.values():
            fh.close()

    specs.sort(key=lambda s: s.shard_index)
    log.info(
        "Built Picard shard manifest: %d shards from %d chroms; "
        "sizes %s",
        len(specs), len(chrom_counts),
        [shard_totals[s.shard_index] for s in specs],
    )
    return specs


def concat_picard_outputs(
    shards: list[PicardShardSpec],
    out_lifted: Path,
    out_rejected: Path,
) -> None:
    """Concat per-shard lifted + rejected VCFs in shard-index order.

    Headers (lines starting with '#') from shard 0 only; body records from all shards.
    Shards are expected to be in shard_index order (build_picard_shard_manifest guarantees this).
    """
    out_lifted.parent.mkdir(parents=True, exist_ok=True)
    for out_path, vcf_attr in [
        (out_lifted, "lifted_vcf"),
        (out_rejected, "rejected_vcf"),
    ]:
        with open(out_path, "w") as out:
            header_written = False
            for shard in shards:
                shard_path: Path = getattr(shard, vcf_attr)
                if not shard_path.exists():
                    raise PileupAadrInternalError(
                        what="concat_picard_outputs",
                        why=f"shard {shard.shard_index} output missing: {shard_path}",
                        fix="Picard should have raised before reaching concat; check stderr",
                    )
                with open(shard_path) as vcf_in:
                    for line in vcf_in:
                        if line.startswith("#"):
                            if not header_written:
                                out.write(line)
                        else:
                            header_written = True
                            out.write(line)


def aggregate_stage1_counters(
    per_shard: list[Stage1LiftCounters],
    input_filters: Stage1InputFilters,
    yield_fail_pct: float,
    yield_warn_pct: float,
) -> Stage1LiftCounters:
    """Field-wise aggregate of per-shard Stage1 counters.

    input_sites_after_filters comes from input_filters.rows_written (the total from
    build_sites_vcf) — not summed from per-shard counters, which each carry the total
    (since lift_aadr_sites is called with the shared input_filters object for each shard).
    wallclock = max (shards run in parallel). liftover_yield_pct recomputed from totals.
    Yield gate evaluated on aggregated totals.
    """
    total_input = input_filters.rows_written
    total_lifted = sum(c.lifted_sites for c in per_shard)
    total_swapped = sum(c.swapped_alleles_count for c in per_shard)
    max_wallclock = max(c.wallclock_seconds for c in per_shard)

    rejected_by_reason: dict[str, int] = {}
    for c in per_shard:
        for reason, count in c.rejected_by_reason.items():
            rejected_by_reason[reason] = rejected_by_reason.get(reason, 0) + count

    yield_pct = 100.0 * total_lifted / total_input if total_input > 0 else 0.0
    yield_warning = yield_pct < yield_warn_pct

    if yield_pct < yield_fail_pct:
        dominant = max(
            rejected_by_reason.items(), key=lambda kv: kv[1], default=("unknown", 0)
        )
        raise LiftoverYieldError(
            what=f"yield {yield_pct:.2f}% (< {yield_fail_pct}% gate)",
            why=(
                f"{total_lifted}/{total_input} sites lifted; dominant rejection: "
                f"{dominant[0]} ({dominant[1]} sites)"
            ),
            fix=(
                "Verify chain file matches BAM build + verify --ref-fasta matches "
                f"BAM @PG. For dominant rejection '{dominant[0]}', see the "
                "troubleshooting guide."
            ),
        )

    if yield_warning:
        dominant_name = max(
            rejected_by_reason.items(), key=lambda kv: kv[1], default=("none", 0)
        )[0]
        log.warning(
            "Liftover yield %.2f%% is below warn threshold %.1f%% (still proceeding); "
            "dominant rejection reason: %s",
            yield_pct, yield_warn_pct, dominant_name,
        )

    log.info(
        "Stage 1 aggregate (%d shards): %d lifted of %d (%.2f%% yield, %d swapped); "
        "max shard wallclock %.1fs",
        len(per_shard), total_lifted, total_input, yield_pct, total_swapped, max_wallclock,
    )

    return Stage1LiftCounters(
        wallclock_seconds=max_wallclock,
        input_sites_after_filters=total_input,
        lifted_sites=total_lifted,
        liftover_yield_pct=round(yield_pct, 4),
        liftover_yield_warning=yield_warning,
        rejected_by_reason=rejected_by_reason,
        swapped_alleles_count=total_swapped,
        input_filters=input_filters,
    )


def lift_aadr_sites_sharded(
    sites_vcf_path: Path,
    chain_path: Path,
    target_fasta_path: Path,
    output_lifted_vcf: Path,
    output_rejected_vcf: Path,
    input_filter_counters: Stage1InputFilters,
    shard_tempdir: Path,
    n_shards: int,
    *,
    picard_mem: str = "3g",
    picard_max_records: int = 100_000,
    yield_fail_pct: float = 70.0,
    yield_warn_pct: float = 95.0,
) -> Stage1LiftCounters:
    """Sharded Picard LiftoverVcf.

    Behavior:
      - n_shards == 1: short-circuit to `lift_aadr_sites` (byte-identical to v0.3).
      - n_shards >  1: partition sites VCF by source chrom (LPT bin-packing),
        run N Picards in parallel, concat lifted+rejected, aggregate counters.

    No-lift fast path interaction: this function is not called when
    `aadr_build == bam_build` — the orchestrator routes around Stage 1 entirely.

    Args:
        sites_vcf_path: sorted sites VCF from build_sites_vcf.
        chain_path: resolved chain file.
        target_fasta_path: TARGET assembly FASTA (hg38; .dict sidecar must exist).
        output_lifted_vcf: merged lifted VCF output path.
        output_rejected_vcf: merged rejected VCF output path.
        input_filter_counters: from build_sites_vcf; stamped onto returned counters unchanged.
        shard_tempdir: parent dir for shard_{00..NN}/ working subdirs.
        n_shards: target Picard parallelism; clamped to chrom count.
        picard_mem: JVM heap per shard (default "3g").
        picard_max_records: --MAX_RECORDS_IN_RAM per shard (default 100_000).
        yield_fail_pct: aggregate yield gate (default 70.0).
        yield_warn_pct: aggregate yield warning (default 95.0).

    Returns:
        Stage1LiftCounters (aggregated; wallclock = max shard time).

    Raises:
        LiftoverYieldError: aggregate yield < yield_fail_pct.
        ToolSubprocessError: any Picard shard exited non-zero (first failure wins).
        PileupAadrInternalError: Picard stderr parse failure.
    """
    if n_shards == 1:
        return lift_aadr_sites(
            sites_vcf_path=sites_vcf_path,
            chain_path=chain_path,
            target_fasta_path=target_fasta_path,
            output_lifted_vcf=output_lifted_vcf,
            output_rejected_vcf=output_rejected_vcf,
            input_filter_counters=input_filter_counters,
            picard_mem=picard_mem,
            picard_max_records=picard_max_records,
            yield_fail_pct=yield_fail_pct,
            yield_warn_pct=yield_warn_pct,
        )

    manifest = build_picard_shard_manifest(sites_vcf_path, shard_tempdir, n_shards)
    actual_shards = len(manifest)

    def _run_shard(spec: PicardShardSpec) -> Stage1LiftCounters:
        # Suppress per-shard yield gates — aggregate after all shards complete.
        # Per-shard denom is wrong (uses total input_filter_counters.rows_written),
        # so skip per-shard gate entirely; aggregate_stage1_counters evaluates correctly.
        return lift_aadr_sites(
            sites_vcf_path=spec.input_vcf,
            chain_path=chain_path,
            target_fasta_path=target_fasta_path,
            output_lifted_vcf=spec.lifted_vcf,
            output_rejected_vcf=spec.rejected_vcf,
            input_filter_counters=input_filter_counters,
            picard_mem=picard_mem,
            picard_max_records=picard_max_records,
            yield_fail_pct=0.0,   # no per-shard gate
            yield_warn_pct=0.0,   # no per-shard warn
        )

    log.info(
        "Stage 1: lifting %d shards in parallel (%d chroms)",
        actual_shards, sum(len(s.source_chroms) for s in manifest),
    )
    per_shard_counters: dict[int, Stage1LiftCounters] = {}
    n_done = 0
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=actual_shards, thread_name_prefix="picard-shard"
    ) as executor:
        futures = {executor.submit(_run_shard, spec): spec for spec in manifest}
        try:
            for future in concurrent.futures.as_completed(futures):
                spec = futures[future]
                per_shard_counters[spec.shard_index] = future.result()
                n_done += 1
                log.info(
                    "Stage 1 shard %d/%d (%s) done: %d sites lifted in %.1fs",
                    n_done, actual_shards,
                    ", ".join(spec.source_chroms),
                    per_shard_counters[spec.shard_index].lifted_sites,
                    per_shard_counters[spec.shard_index].wallclock_seconds,
                )
        except Exception:
            executor.shutdown(wait=True, cancel_futures=True)
            raise

    concat_picard_outputs(manifest, output_lifted_vcf, output_rejected_vcf)

    # Combine per-shard stderrs to the canonical location for debugging
    combined_stderr = output_lifted_vcf.parent / "picard.stderr"
    with open(combined_stderr, "w") as out_stderr:
        for shard in manifest:
            if shard.stderr_path.exists():
                out_stderr.write(f"# --- shard {shard.shard_index} ({', '.join(shard.source_chroms)}) ---\n")
                out_stderr.write(shard.stderr_path.read_text())

    ordered_counters = [per_shard_counters[s.shard_index] for s in manifest]
    return aggregate_stage1_counters(
        per_shard=ordered_counters,
        input_filters=input_filter_counters,
        yield_fail_pct=yield_fail_pct,
        yield_warn_pct=yield_warn_pct,
    )


def lift_and_transform_sharded(
    sites_vcf_path: Path,
    chain_path: Path,
    target_fasta_path: Path,
    output_lifted_vcf: Path,
    output_rejected_vcf: Path,
    output_snp_path: Path,
    output_bed_path: Path,
    input_filter_counters: Stage1InputFilters,
    shard_tempdir: Path,
    n_shards: int,
    *,
    alt_contig_filter: bool = True,
    picard_mem: str = "3g",
    picard_max_records: int = 100_000,
    yield_fail_pct: float = 70.0,
    yield_warn_pct: float = 95.0,
) -> tuple[Stage1LiftCounters, Stage2TransformCounters]:
    """Streaming Stage 1 + Stage 2: transform each shard as soon as Picard finishes it.

    n_shards == 1: sequential lift then transform (byte-identical to calling them separately).
    n_shards >  1: Picard shards run in parallel; as each shard completes its lifted VCF
        is immediately submitted for Stage 2 transform while remaining Picard shards run.
        Eliminates transform serialization latency (~transform_time * (1 - 1/n_shards)).

    The concatenated lifted VCF (`output_lifted_vcf`) is still written for Stage 4's
    swap-lookup build. Per-shard .snp/.bed fragments are concatenated into
    `output_snp_path` / `output_bed_path` in shard_index order; Stage 3's
    `build_shard_manifest` re-groups by chromosome, so global sort order is not required.
    """
    from . import transform as _transform  # local import avoids module-level cycle risk

    if n_shards == 1:
        s1 = lift_aadr_sites(
            sites_vcf_path=sites_vcf_path,
            chain_path=chain_path,
            target_fasta_path=target_fasta_path,
            output_lifted_vcf=output_lifted_vcf,
            output_rejected_vcf=output_rejected_vcf,
            input_filter_counters=input_filter_counters,
            picard_mem=picard_mem,
            picard_max_records=picard_max_records,
            yield_fail_pct=yield_fail_pct,
            yield_warn_pct=yield_warn_pct,
        )
        s2 = _transform.build_pileupcaller_snp_and_bed(
            lifted_vcf_path=output_lifted_vcf,
            output_snp_path=output_snp_path,
            output_bed_path=output_bed_path,
            alt_contig_filter=alt_contig_filter,
        )
        return s1, s2

    manifest = build_picard_shard_manifest(sites_vcf_path, shard_tempdir, n_shards)
    actual_shards = len(manifest)
    log.info(
        "Stage 1+2: lifting %d shards in parallel (%d chroms)",
        actual_shards, sum(len(s.source_chroms) for s in manifest),
    )

    def _run_shard(spec: PicardShardSpec) -> Stage1LiftCounters:
        return lift_aadr_sites(
            sites_vcf_path=spec.input_vcf,
            chain_path=chain_path,
            target_fasta_path=target_fasta_path,
            output_lifted_vcf=spec.lifted_vcf,
            output_rejected_vcf=spec.rejected_vcf,
            input_filter_counters=input_filter_counters,
            picard_mem=picard_mem,
            picard_max_records=picard_max_records,
            yield_fail_pct=0.0,
            yield_warn_pct=0.0,
        )

    def _run_transform(spec: PicardShardSpec) -> Stage2TransformCounters:
        snp_frag = spec.lifted_vcf.parent / "transform.snp"
        bed_frag = spec.lifted_vcf.parent / "transform.bed"
        return _transform.build_pileupcaller_snp_and_bed(
            lifted_vcf_path=spec.lifted_vcf,
            output_snp_path=snp_frag,
            output_bed_path=bed_frag,
            alt_contig_filter=alt_contig_filter,
        )

    per_shard_s1: dict[int, Stage1LiftCounters] = {}
    per_shard_s2: dict[int, Stage2TransformCounters] = {}
    transform_futures: dict[concurrent.futures.Future[Stage2TransformCounters], PicardShardSpec] = {}
    n_lift_done = 0
    n_transform_done = 0

    with concurrent.futures.ThreadPoolExecutor(
        max_workers=actual_shards, thread_name_prefix="picard-shard"
    ) as executor:
        picard_futures = {executor.submit(_run_shard, spec): spec for spec in manifest}
        try:
            for picard_future in concurrent.futures.as_completed(picard_futures):
                spec = picard_futures[picard_future]
                per_shard_s1[spec.shard_index] = picard_future.result()
                n_lift_done += 1
                log.info(
                    "Stage 1 shard %d/%d (%s) done: %d sites lifted in %.1fs — queuing Stage 2",
                    n_lift_done, actual_shards, ", ".join(spec.source_chroms),
                    per_shard_s1[spec.shard_index].lifted_sites,
                    per_shard_s1[spec.shard_index].wallclock_seconds,
                )
                transform_futures[executor.submit(_run_transform, spec)] = spec
        except Exception:
            executor.shutdown(wait=True, cancel_futures=True)
            raise

        try:
            for t_future in concurrent.futures.as_completed(transform_futures):
                spec = transform_futures[t_future]
                per_shard_s2[spec.shard_index] = t_future.result()
                n_transform_done += 1
                log.info(
                    "Stage 2 shard %d/%d (%s) done: %d sites written in %.1fs",
                    n_transform_done, actual_shards, ", ".join(spec.source_chroms),
                    per_shard_s2[spec.shard_index].output_sites,
                    per_shard_s2[spec.shard_index].wallclock_seconds,
                )
        except Exception:
            executor.shutdown(wait=True, cancel_futures=True)
            raise

    concat_picard_outputs(manifest, output_lifted_vcf, output_rejected_vcf)

    combined_stderr = output_lifted_vcf.parent / "picard.stderr"
    with open(combined_stderr, "w") as out_stderr:
        for shard in manifest:
            if shard.stderr_path.exists():
                out_stderr.write(
                    f"# --- shard {shard.shard_index} "
                    f"({', '.join(shard.source_chroms)}) ---\n"
                )
                out_stderr.write(shard.stderr_path.read_text())

    output_snp_path.parent.mkdir(parents=True, exist_ok=True)
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_snp_path, "w") as snp_out, open(output_bed_path, "w") as bed_out:
        for shard in manifest:
            snp_frag = shard.lifted_vcf.parent / "transform.snp"
            bed_frag = shard.lifted_vcf.parent / "transform.bed"
            if snp_frag.exists():
                snp_out.write(snp_frag.read_text())
            if bed_frag.exists():
                bed_out.write(bed_frag.read_text())

    ordered_s1 = [per_shard_s1[s.shard_index] for s in manifest]
    s1 = aggregate_stage1_counters(
        per_shard=ordered_s1,
        input_filters=input_filter_counters,
        yield_fail_pct=yield_fail_pct,
        yield_warn_pct=yield_warn_pct,
    )

    ordered_s2 = [per_shard_s2[s.shard_index] for s in manifest]
    s2 = Stage2TransformCounters(
        wallclock_seconds=max(c.wallclock_seconds for c in ordered_s2),
        alt_contig_drops=sum(c.alt_contig_drops for c in ordered_s2),
        output_sites=sum(c.output_sites for c in ordered_s2),
    )

    log.info(
        "Stage 1+2 streaming complete (%d shards): %d sites; "
        "Stage 2 max shard wallclock %.2fs",
        actual_shards, s2.output_sites, s2.wallclock_seconds,
    )
    return s1, s2


def parse_picard_stderr(stderr_text: str) -> dict[str, int]:
    """Extract structured counters from Picard LiftoverVcf stderr.

    Returns:
        {
            "processed_count": int,
            "failed_count": int,
            "mismatched_ref_count": int,
            "lifted_count": int,             # = processed - failed
            "swapped_count": int,             # SwappedAlleles INFO marker count
        }

    Raises:
        PileupAadrInternalError: a REQUIRED pattern fails to match (Picard
            stderr format change). Optional patterns (currently just `swapped`)
            use their default value silently.
    """
    extracted: dict[str, str] = {}
    for key, pattern in _PICARD_REQUIRED_PATTERNS.items():
        m = pattern.search(stderr_text)
        if m is None:
            raise PileupAadrInternalError(
                what="parse_picard_stderr",
                why=f"failed to match required pattern {key!r}: {pattern.pattern!r}",
                fix=(
                    "Picard stderr format may have changed. Verify Picard version "
                    f"matches PICARD_SPEC.tested_against ({PICARD_SPEC.tested_against}); "
                    "if format truly changed, file a bug report with the stderr output attached"
                ),
            )
        extracted[key] = m.group(1)

    for key, (pattern, default) in _PICARD_OPTIONAL_PATTERNS.items():
        m = pattern.search(stderr_text)
        if m is None:
            log.debug("Picard stderr lacks optional pattern %r; defaulting to %r", key, default)
            extracted[key] = default
        else:
            extracted[key] = m.group(1)

    processed = int(extracted["processed"])
    failed = int(extracted["failed"])
    return {
        "processed_count": processed,
        "failed_count": failed,
        "mismatched_ref_count": int(extracted["mismatched_ref"]),
        "lifted_count": processed - failed,
        "swapped_count": int(extracted["swapped"]),
    }


def parse_rejected_vcf(rejected_vcf_path: Path) -> dict[str, int]:
    """Parse Picard's `_rejected.vcf` to count by FILTER reason.

    Empty rejected.vcf is fine; returns all-zero dict.

    Returns:
        {
            "NoTarget": N,
            "MismatchedRefAllele": N,
            "IndelStraddlesMultipleIntervals": N,
            "SwappedAlleles": N,    # should always be 0 (recovered, not rejected)
            "other": N,
        }
    """
    counts: dict[str, int] = dict.fromkeys(_REJECT_FILTER_BUCKETS, 0)
    counts["other"] = 0
    if not rejected_vcf_path.exists():
        return counts

    with pysam.VariantFile(str(rejected_vcf_path)) as vcf:
        for record in vcf:
            filters = list(record.filter.keys())
            if not filters or filters == ["PASS"]:
                # Shouldn't happen in _rejected.vcf but defensive
                counts["other"] += 1
                continue
            primary = filters[0]
            if primary in counts:
                counts[primary] += 1
            else:
                counts["other"] += 1
                log.debug("Unknown Picard rejection FILTER %r", primary)
    return counts


__all__ = [
    "PicardShardSpec",
    "aggregate_stage1_counters",
    "build_picard_shard_manifest",
    "chain_file_path",
    "concat_picard_outputs",
    "get_bundled_chain_path",
    "lift_aadr_sites",
    "lift_aadr_sites_sharded",   # public: Stage 1 only (no transform); orchestrator uses lift_and_transform_sharded
    "lift_and_transform_sharded",
    "parse_picard_stderr",
    "parse_rejected_vcf",
    "resolve_chain_for_extract",
]
