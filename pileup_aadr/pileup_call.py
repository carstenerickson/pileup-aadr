"""Stage 3: pipe samtools mpileup output into pileupCaller --randomDiploid.

Longest wallclock stage (~30-40 min on a 33x WGS at 1240k). Two subprocesses
connected by an OS pipe via `ToolWrapper.pipe`. SIGPIPE handling: if the
downstream (pileupCaller) dies first, samtools exits 141 (128 + SIGPIPE);
we tolerate upstream 141 IFF the downstream exited cleanly.

Module name choice (v2.2 M12): named `pileup_call` rather than `call` to avoid
the unqualified `call` import alias that would clash with builtin terminology.
"""
from __future__ import annotations

import concurrent.futures
import logging
import re
import shutil
import time
from pathlib import Path
from typing import Final

from .counters import PileupCallerSummary, Stage3CallCounters, Stage3ShardCounters
from .errors import PileupAadrInternalError, ToolSubprocessError
from .shard import ShardSpec, build_shard_manifest, merge_shard_eigenstrat
from .tool_wrapper import PILEUPCALLER_SPEC, SAMTOOLS_SPEC, ToolWrapper

log = logging.getLogger(__name__)

# pileupCaller 1.6.0.0 stderr summary block (verified empirically v2.1).
# Format:
#   # Summary Statistics per sample
#   # ...
#   SampleName\tTotalSites\tNonMissingCalls\tavgRawReads\tavgDamageCleanedReads\tavgSampledFrom
#   <name>\t<int>\t<int>\t<float>\t<float>\t<float>
_PILEUPCALLER_HEADER_RE: Final[re.Pattern[str]] = re.compile(
    r"^SampleName\s+TotalSites\s+NonMissingCalls"
    r"\s+avgRawReads\s+avgDamageCleanedReads\s+avgSampledFrom",
    re.MULTILINE,
)
# Scientific notation is emitted at low coverage (e.g., "3.68e-2") so floats
# need the full pattern, not just `[\d.]+`. The integer columns stay strict.
_NUM_RE: Final[str] = r"\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"
_PILEUPCALLER_DATA_RE: Final[re.Pattern[str]] = re.compile(
    rf"^(\S+)\s+(\d+)\s+(\d+)\s+({_NUM_RE})\s+({_NUM_RE})\s+({_NUM_RE})\s*$",
    re.MULTILINE,
)

# SIGPIPE-induced upstream exit code (128 + 13). Tolerated only when downstream
# exited cleanly (pipe-shutdown ordering).
_SIGPIPE_EXIT: Final[int] = 141


def run_pileup_call(
    bam_path: Path,
    snp_path: Path,
    bed_path: Path,
    target_fasta_path: Path,
    output_prefix: Path,
    sample_name: str,
    pop_name: str,
    *,
    seed: int = 42,
    min_mapq: int = 30,
    min_baseq: int = 30,
    no_baq: bool = False,
    region: str | None = None,
) -> Stage3CallCounters:
    """Pipe `samtools mpileup` output into `pileupCaller --randomDiploid`.

    Args:
        bam_path: aligned BAM/CRAM (BAM index .bai/.crai must exist alongside).
        snp_path: pileupCaller .snp from Stage 2.
        bed_path: mpileup BED from Stage 2.
        target_fasta_path: TARGET assembly FASTA (must match BAM build).
        output_prefix: pileupCaller writes `<prefix>.{geno,snp,ind}` here.
        sample_name: --sampleNames value; becomes IID in output .ind.
        pop_name: --samplePopName value; becomes POP in output .ind.
        seed: pileupCaller --seed (default 42 for reproducibility).
        min_mapq: mpileup -q (default 30).
        min_baseq: mpileup -Q (default 30).
        no_baq: disable mpileup -B (default: -B on; pileup-aadr is for
            modern-WGS-→-AADR-projection so -B helps with REF bias).
        region: samtools region string (e.g. "chr1") for per-shard calls.
            None means no region restriction (full-BAM single-process path).

    Returns:
        Stage3CallCounters with parsed pileupCaller stderr summary + wallclock.
        per_shard is always empty — call-level counters only.

    Raises:
        ToolSubprocessError: samtools or pileupCaller exited non-zero (modulo
            tolerated SIGPIPE on samtools when pileupCaller succeeded).
        PileupAadrInternalError: pileupCaller stderr unparseable (format change).
    """
    samtools_args = ["mpileup"]
    if not no_baq:
        samtools_args.append("-B")
    samtools_args.extend([
        f"-q{min_mapq}",
        f"-Q{min_baseq}",
        "-R",  # ignore-RG (one BAM = one sample); per pileupCaller's recommended call
        "-f", str(target_fasta_path),
        "-l", str(bed_path),
    ])
    samtools_args.append(str(bam_path))
    if region is not None:
        samtools_args.append(region)

    pileupcaller_args = [
        "--randomDiploid",
        "--seed", str(seed),
        "-f", str(snp_path),
        "--sampleNames", sample_name,
        "--samplePopName", pop_name,
        "-e", str(output_prefix),
    ]

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    mpileup_stderr = output_prefix.parent / "mpileup.stderr"
    pileupcaller_stderr = output_prefix.parent / "pileupcaller.stderr"

    samtools_wrapper = ToolWrapper(SAMTOOLS_SPEC)
    pileupcaller_wrapper = ToolWrapper(PILEUPCALLER_SPEC)

    t0 = time.perf_counter()
    upstream_result, downstream_result = samtools_wrapper.pipe(
        downstream=pileupcaller_wrapper,
        upstream_args=samtools_args,
        downstream_args=pileupcaller_args,
        upstream_stderr_to=mpileup_stderr,
        downstream_stderr_to=pileupcaller_stderr,
    )
    wallclock = time.perf_counter() - t0

    if downstream_result.exit_code != 0:
        raise ToolSubprocessError(
            what=f"pileupCaller (stderr at {pileupcaller_stderr})",
            why=(
                f"exit code {downstream_result.exit_code}; tail of stderr:\n"
                f"{_tail(pileupcaller_stderr, 20)}"
            ),
            fix=(
                "Check stderr for pileupCaller's diagnostic; common causes: "
                "alt-contig in .snp file (Stage 2 alt_contig_filter should have "
                "caught), chromosome-name mismatch between .snp and BAM"
            ),
        )
    if upstream_result.exit_code not in (0, _SIGPIPE_EXIT):
        raise ToolSubprocessError(
            what=f"samtools mpileup (stderr at {mpileup_stderr})",
            why=(
                f"exit code {upstream_result.exit_code}; tail of stderr:\n"
                f"{_tail(mpileup_stderr, 20)}"
            ),
            fix=(
                "Check stderr for samtools' diagnostic; common causes: "
                "missing BAM index, FASTA/.fai mismatch, BED format error"
            ),
        )

    summary = parse_pileupcaller_stderr(pileupcaller_stderr.read_text())

    log.debug(
        "Shard %s complete: pileupCaller called %d / %d sites "
        "(%.1fx avg coverage); wallclock %.1fs",
        region or "all",
        summary.non_missing_calls, summary.total_sites,
        summary.avg_raw_reads, wallclock,
    )

    return Stage3CallCounters(
        wallclock_seconds=wallclock,
        pileupcaller_summary=summary,
    )


def run_pileup_call_shards(
    bam_path: Path,
    sites_snp_path: Path,
    sites_bed_path: Path,
    target_fasta_path: Path,
    output_prefix: Path,
    sample_name: str,
    pop_name: str,
    shard_dir: Path,
    master_seed: int,
    *,
    threads: int = 1,
    min_mapq: int = 30,
    min_baseq: int = 30,
    no_baq: bool = False,
) -> Stage3CallCounters:
    """Run Stage 3 across N parallel per-chromosome shards.

    When threads=1, short-circuits to a single run_pileup_call invocation
    (no fan-out, no shard manifest, no merge) — byte-identical to v0.2 output.

    When threads>1, partitions the sites into per-chromosome shards, runs up to
    `threads` shards concurrently via ThreadPoolExecutor, then concatenates the
    per-shard EIGENSTRAT triplets in CHROM_ORDER before returning.

    Args:
        bam_path: aligned BAM/CRAM.
        sites_snp_path: Stage 2 output .snp (feeds shard manifest + per-shard call).
        sites_bed_path: Stage 2 output BED (same).
        target_fasta_path: TARGET assembly FASTA.
        output_prefix: merged EIGENSTRAT prefix (<prefix>.{geno,snp,ind}).
        sample_name: --sampleNames value.
        pop_name: --samplePopName value.
        shard_dir: directory for per-chromosome shard working files.
        master_seed: pileupCaller master seed; per-shard seeds derived via
            master_seed * 1009 + shard_index.
        threads: parallelism width (1 = no fan-out; default 1).
        min_mapq: mpileup -q.
        min_baseq: mpileup -Q.
        no_baq: disable mpileup -B.

    Returns:
        Stage3CallCounters. per_shard is populated when threads>1; empty when
        threads==1 (no fan-out decomposition).

    Raises:
        ToolSubprocessError, PileupAadrInternalError: propagated from the first
            failing shard (remaining shards cancelled on first failure).
    """
    if threads == 1:
        counters = run_pileup_call(
            bam_path=bam_path,
            snp_path=sites_snp_path,
            bed_path=sites_bed_path,
            target_fasta_path=target_fasta_path,
            output_prefix=output_prefix,
            sample_name=sample_name,
            pop_name=pop_name,
            seed=master_seed,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            no_baq=no_baq,
        )
        log.info(
            "Stage 3 complete (single-process): pileupCaller called %d / %d sites "
            "(%.1fx avg coverage); wallclock %.1fs",
            counters.pileupcaller_summary.non_missing_calls,
            counters.pileupcaller_summary.total_sites,
            counters.pileupcaller_summary.avg_raw_reads,
            counters.wallclock_seconds,
        )
        return counters

    manifest = build_shard_manifest(
        sites_snp_path=sites_snp_path,
        sites_bed_path=sites_bed_path,
        shard_dir=shard_dir,
        master_seed=master_seed,
    )

    wall_t0 = time.perf_counter()
    shard_results: dict[int, Stage3CallCounters] = {}

    def _run_shard(spec: ShardSpec) -> tuple[ShardSpec, Stage3CallCounters]:
        return spec, run_pileup_call(
            bam_path=bam_path,
            snp_path=spec.snp_path,
            bed_path=spec.bed_path,
            target_fasta_path=target_fasta_path,
            output_prefix=spec.output_prefix,
            sample_name=sample_name,
            pop_name=pop_name,
            seed=spec.seed,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            no_baq=no_baq,
            region=spec.chromosome,
        )

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(_run_shard, spec): spec for spec in manifest}
        try:
            for future in concurrent.futures.as_completed(futures):
                spec, shard_counters = future.result()
                shard_results[spec.shard_index] = shard_counters
                log.debug(
                    "Shard %s done: %d / %d sites; wallclock %.1fs",
                    spec.chromosome,
                    shard_counters.pileupcaller_summary.non_missing_calls,
                    shard_counters.pileupcaller_summary.total_sites,
                    shard_counters.wallclock_seconds,
                )
        except Exception:
            for f in futures:
                f.cancel()
            raise

    merge_shard_eigenstrat(manifest, output_prefix)

    for shard in manifest:
        shutil.rmtree(shard.bed_path.parent, ignore_errors=True)
        log.debug("Cleaned shard dir: %s", shard.bed_path.parent)

    wall_total = time.perf_counter() - wall_t0

    # Aggregate per-shard summaries (call-count-weighted for avg fields)
    total_sites = 0
    total_non_missing = 0
    w_raw = w_dc = w_sf = 0.0
    for s in manifest:
        pc = shard_results[s.shard_index].pileupcaller_summary
        total_sites += pc.total_sites
        total_non_missing += pc.non_missing_calls
        w_raw += pc.avg_raw_reads * pc.total_sites
        w_dc += pc.avg_damage_cleaned_reads * pc.total_sites
        w_sf += pc.avg_sampled_from * pc.total_sites
    if total_sites > 0:
        avg_raw, avg_dc, avg_sf = w_raw / total_sites, w_dc / total_sites, w_sf / total_sites
    else:
        avg_raw = avg_dc = avg_sf = 0.0

    per_shard = [
        Stage3ShardCounters(
            shard_index=spec.shard_index,
            chromosome=spec.chromosome,
            pileupcaller_summary=shard_results[spec.shard_index].pileupcaller_summary,
            wallclock_seconds=shard_results[spec.shard_index].wallclock_seconds,
        )
        for spec in manifest
    ]

    aggregate_summary = PileupCallerSummary(
        total_sites=total_sites,
        non_missing_calls=total_non_missing,
        avg_raw_reads=avg_raw,
        avg_damage_cleaned_reads=avg_dc,
        avg_sampled_from=avg_sf,
    )

    log.info(
        "Stage 3 complete (%d shards, %d threads): pileupCaller called %d / %d sites "
        "(%.1fx avg coverage); wallclock %.1fs",
        len(manifest), threads,
        aggregate_summary.non_missing_calls, aggregate_summary.total_sites,
        aggregate_summary.avg_raw_reads, wall_total,
    )

    return Stage3CallCounters(
        wallclock_seconds=wall_total,
        pileupcaller_summary=aggregate_summary,
        per_shard=per_shard,
    )


def parse_pileupcaller_stderr(stderr_text: str) -> PileupCallerSummary:
    """Parse pileupCaller's structured stderr summary stats TSV.

    Format (single sample — pileup-aadr is single-sample by design):

        # Summary Statistics per sample
        # ...
        SampleName  TotalSites  NonMissingCalls  avgRawReads  avgDamageCleanedReads  avgSampledFrom
        Carsten     1131600     1107213          21.4         21.4                   21.4

    Raises:
        PileupAadrInternalError: stderr lacks the expected header or data line
            (pileupCaller version-skew → fail closed).
    """
    header_match = _PILEUPCALLER_HEADER_RE.search(stderr_text)
    if header_match is None:
        raise PileupAadrInternalError(
            what="parse_pileupcaller_stderr",
            why="expected summary-stats header line not found in stderr",
            fix=(
                "pileupCaller stderr format may have changed. Verify pileupCaller "
                f"version matches PILEUPCALLER_SPEC.tested_against "
                f"({PILEUPCALLER_SPEC.tested_against}); if format truly changed, "
                "file a bug report with the stderr output attached"
            ),
        )
    data_match = _PILEUPCALLER_DATA_RE.search(stderr_text, pos=header_match.end())
    if data_match is None:
        raise PileupAadrInternalError(
            what="parse_pileupcaller_stderr",
            why="expected data line (sample row) not found after header",
            fix=(
                "pileupCaller stderr format may have changed; file a bug report "
                "with stderr attached"
            ),
        )
    return PileupCallerSummary(
        total_sites=int(data_match.group(2)),
        non_missing_calls=int(data_match.group(3)),
        avg_raw_reads=float(data_match.group(4)),
        avg_damage_cleaned_reads=float(data_match.group(5)),
        avg_sampled_from=float(data_match.group(6)),
    )


def _tail(path: Path, n: int) -> str:
    if not path.exists():
        return "(stderr file missing)"
    return "\n".join(path.read_text().splitlines()[-n:])


__all__ = ["parse_pileupcaller_stderr", "run_pileup_call", "run_pileup_call_shards"]
