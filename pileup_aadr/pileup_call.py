"""Stage 3: pipe samtools mpileup output into pileupCaller --randomDiploid.

Longest wallclock stage (~30-40 min on a 33x WGS at 1240k). Two subprocesses
connected by an OS pipe via `ToolWrapper.pipe`. SIGPIPE handling: if the
downstream (pileupCaller) dies first, samtools exits 141 (128 + SIGPIPE);
we tolerate upstream 141 IFF the downstream exited cleanly.

Module name choice (v2.2 M12): named `pileup_call` rather than `call` to avoid
the unqualified `call` import alias that would clash with builtin terminology.
"""
from __future__ import annotations

import logging
import re
import time
from pathlib import Path
from typing import Final

from .counters import PileupCallerSummary, Stage3CallCounters
from .errors import PileupAadrInternalError, ToolSubprocessError
from .tool_wrapper import PILEUPCALLER_SPEC, SAMTOOLS_SPEC, ToolWrapper

log = logging.getLogger(__name__)

# mpileup is BAM-seek-bound on large BAMs; >4 threads gives diminishing returns
# (verified empirically v2.1 — wallclock plateau at threads=4 for 67 GB BAM).
_THREAD_CAP: Final[int] = 4

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
    threads: int = 1,
    min_mapq: int = 30,
    min_baseq: int = 30,
    no_baq: bool = False,
    no_thread_cap: bool = False,
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
        threads: samtools mpileup -@ N (auto-capped to 4 unless `no_thread_cap`).
        min_mapq: mpileup -q (default 30).
        min_baseq: mpileup -Q (default 30).
        no_baq: disable mpileup -B (default: -B on; pileup-aadr is for
            modern-WGS-→-AADR-projection so -B helps with REF bias).
        no_thread_cap: skip the threads cap (default False; cap is on).

    Returns:
        Stage3CallCounters with parsed pileupCaller stderr summary + wallclock.

    Raises:
        ToolSubprocessError: samtools or pileupCaller exited non-zero (modulo
            tolerated SIGPIPE on samtools when pileupCaller succeeded).
        PileupAadrInternalError: pileupCaller stderr unparseable (format change).
    """
    effective_threads = threads if no_thread_cap else min(threads, _THREAD_CAP)
    if effective_threads != threads:
        log.info(
            "Capping mpileup threads %d -> %d (mpileup is BAM-seek-bound on large "
            "BAMs; pass --no-thread-cap to override)",
            threads, effective_threads,
        )

    # samtools mpileup args
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
    if effective_threads > 1:
        samtools_args.extend(["-@", str(effective_threads)])
    samtools_args.append(str(bam_path))

    # pileupCaller args
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

    # Check downstream first — its output is what matters
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
    # Upstream 141 (SIGPIPE) tolerated only if downstream succeeded (verified above).
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

    log.info(
        "Stage 3 complete: pileupCaller called %d / %d sites "
        "(%.1fx avg coverage); wallclock %.1fs",
        summary.non_missing_calls, summary.total_sites,
        summary.avg_raw_reads, wallclock,
    )

    return Stage3CallCounters(
        wallclock_seconds=wallclock,
        pileupcaller_summary=summary,
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
    if _PILEUPCALLER_HEADER_RE.search(stderr_text) is None:
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
    header_match = _PILEUPCALLER_HEADER_RE.search(stderr_text)
    assert header_match is not None
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


__all__ = ["parse_pileupcaller_stderr", "run_pileup_call"]
