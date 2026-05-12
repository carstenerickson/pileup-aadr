"""`validate` subcommand implementation — 10 pre-flight checks per LLD §16.

Day-1 scope: file-side checks fully implemented (BAM index, BAM build, AADR parse, AADR
build, output-prefix collisions, ref-FASTA findability). Binary-version checks are stubbed
to PASS-with-skip-reason because `tool_wrapper.py` lands on Day 2; the stubs are tagged
with TODO comments naming the Day-2 implementation gap.
"""
from __future__ import annotations

import json
import logging
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from .errors import PileupAadrError
from .format_detect import (
    BuildOverride,
    detect_aadr_build,
    detect_bam_build,
    detect_bam_format,
    parse_aadr_snp,
)

log = logging.getLogger(__name__)


@dataclass
class CheckResult:
    """One pre-flight check result."""

    name: str
    status: Literal["PASS", "WARN", "FAIL", "SKIP"]
    detail: str = ""


def run_validate(
    *,
    bam: Path,
    aadr_snp: Path,
    chain_path: Path | None,
    ref_fasta: Path | None,
    bam_build_override: BuildOverride,
    aadr_build_override: BuildOverride,
    json_output: bool,
    output_prefix: Path | None,
) -> int:
    """Run the 10 pre-flight checks. Returns exit code (0 if all PASS, 1 if any FAIL).

    Checks per LLD §16:
        Tool dependencies:
            1. samtools binary present + version >= 1.16
            2. picard binary present + version >= 3.0 (skipped if no-lift fast path)
            3. pileupCaller binary present + version >= 1.6.0
            4. java present + version >= 11 (skipped if picard not needed)
        Input files:
            5. BAM/CRAM readable + index present
            6. BAM build detectable
            7. AADR .snp parseable + no duplicate rsIDs
            8. AADR build detectable
            9. Target FASTA findable + .dict OK + chrom-name match BAM @SQ
        Output:
            10. Output prefix not lock-held + parent directory writable (skipped if -o not given)
    """
    results: list[CheckResult] = []

    # --- AADR-side detection first (informs whether picard/java are needed) ---
    aadr_df = None
    aadr_build = None
    try:
        aadr_df = parse_aadr_snp(aadr_snp)
        results.append(_aadr_parse_pass(aadr_df))
    except PileupAadrError as e:
        results.append(CheckResult("AADR .snp parseable", "FAIL", e.why))

    if aadr_df is not None:
        try:
            aadr_build = detect_aadr_build(aadr_df, override=aadr_build_override)
            results.append(
                CheckResult("AADR build detectable", "PASS", f"detected: {aadr_build}")
            )
        except PileupAadrError as e:
            results.append(CheckResult("AADR build detectable", "FAIL", e.why))

    # --- BAM-side detection ---
    bam_build = None
    try:
        bam_format = detect_bam_format(bam)
        results.append(_bam_index_check(bam, bam_format))
        bam_build = detect_bam_build(bam, override=bam_build_override)
        results.append(
            CheckResult("BAM build detectable", "PASS", f"{bam_format} {bam_build}")
        )
    except PileupAadrError as e:
        results.append(CheckResult("BAM readable + indexed", "FAIL", e.why))

    # --- determine no-lift fast path ---
    no_lift = bam_build is not None and aadr_build is not None and bam_build == aadr_build

    # --- Tool dependencies ---
    # Day-1 scope note: real version checks land on Day 2 (tool_wrapper.py).
    # For now, we just confirm the binary is on PATH and report SKIP-with-reason if not.
    results.append(_check_binary_present("samtools", required=True))
    if no_lift:
        results.append(CheckResult("picard binary", "SKIP", "no-lift fast path"))
        results.append(CheckResult("java binary", "SKIP", "no-lift fast path"))
    else:
        results.append(_check_binary_present("picard", required=True, also_check=("picard.jar",)))
        results.append(_check_binary_present("java", required=True))
    results.append(_check_binary_present("pileupCaller", required=True))

    # --- Chain ---
    if chain_path is not None and not chain_path.exists():
        results.append(
            CheckResult("chain file", "FAIL", f"--chain {chain_path} not found")
        )
    else:
        results.append(
            CheckResult(
                "chain file",
                "PASS",
                "user-supplied" if chain_path is not None else "bundled (default)",
            )
        )

    # --- Target FASTA ---
    results.append(_check_ref_fasta(ref_fasta, bam, bam_build))

    # --- Output prefix ---
    if output_prefix is None:
        results.append(CheckResult("output prefix", "SKIP", "no -o passed"))
    else:
        results.append(_check_output_prefix(output_prefix))

    return _emit_results(results, json_output=json_output)


def _aadr_parse_pass(aadr_df: object) -> CheckResult:
    """parse_aadr_snp succeeded — render a PASS with summary."""
    import pandas as pd

    if not isinstance(aadr_df, pd.DataFrame):
        return CheckResult("AADR .snp parseable", "FAIL", "internal error: non-DataFrame")
    return CheckResult(
        "AADR .snp parseable",
        "PASS",
        f"{len(aadr_df)} unique rsIDs, {aadr_df['chrom_int'].nunique()} chromosomes",
    )


def _bam_index_check(bam: Path, bam_format: str) -> CheckResult:
    """Check that the BAM has a sibling .bai (or .crai for CRAM)."""
    if bam_format == "BAM":
        bai_candidates = (Path(f"{bam}.bai"), bam.with_suffix(".bai"))
    else:
        bai_candidates = (Path(f"{bam}.crai"), bam.with_suffix(".crai"))
    if not any(p.exists() for p in bai_candidates):
        return CheckResult(
            "BAM index", "FAIL", f"no index ({'/'.join(p.suffix for p in bai_candidates)}) found"
        )
    return CheckResult("BAM readable + indexed", "PASS", f"{bam_format}, index present")


def _check_binary_present(
    binary: str, *, required: bool, also_check: tuple[str, ...] = ()
) -> CheckResult:
    """Check that a binary is on PATH or available at one of the also_check paths.

    Day-1 stub: just on-PATH presence; version check lands on Day 2 via tool_wrapper.
    """
    found_path = shutil.which(binary)
    if found_path is not None:
        return CheckResult(
            f"{binary} binary",
            "PASS",
            f"{found_path} (Day-2 TODO: version-check via tool_wrapper)",
        )
    # Try fallback locations (e.g., picard.jar in conda paths)
    for candidate in also_check:
        if shutil.which(candidate) is not None:
            return CheckResult(
                f"{binary} binary",
                "PASS",
                f"{candidate} found (Day-2 TODO: version-check)",
            )
    return CheckResult(
        f"{binary} binary",
        "FAIL" if required else "WARN",
        f"not found on PATH; install via `conda install -c bioconda {binary.lower()}`",
    )


def _check_ref_fasta(
    cli_ref: Path | None, bam: Path | None, bam_build: str | None
) -> CheckResult:
    """Resolve target FASTA per the LLD §15 resolution order; report what was found."""
    if cli_ref is not None:
        if cli_ref.exists():
            return CheckResult("ref FASTA findable", "PASS", f"--ref-fasta: {cli_ref}")
        return CheckResult(
            "ref FASTA findable", "FAIL", f"--ref-fasta {cli_ref} not found"
        )
    env_dir = os.environ.get("PILEUP_AADR_REF_DIR")
    if env_dir and bam_build is not None:
        candidate = Path(env_dir) / f"{bam_build}.fa"
        if candidate.exists():
            return CheckResult(
                "ref FASTA findable", "PASS", f"$PILEUP_AADR_REF_DIR: {candidate}"
            )
    # BAM @PG fallback — full check requires extract_orch's _extract_ref_from_bam_pg
    # which lands on Day 3. For Day 1 just report SKIP-with-reason.
    return CheckResult(
        "ref FASTA findable",
        "SKIP",
        "Day-2 TODO: BAM @PG fallback lookup",
    )


def _check_output_prefix(prefix: Path) -> CheckResult:
    """Check parent dir writable + no existing output files + no active lock holder.

    Per LLD v2.1 H10 fix: do NOT acquire the lock here (gives false confidence about
    future acquireability); just check for an existing holder via the sidecar PID file.
    """
    if not prefix.parent.exists():
        return CheckResult(
            "output prefix writable",
            "FAIL",
            f"parent dir {prefix.parent} missing",
        )
    if not os.access(prefix.parent, os.W_OK):
        return CheckResult(
            "output prefix writable",
            "FAIL",
            f"parent dir {prefix.parent} not writable",
        )

    # Holder sidecar check (concurrency.py lands on Day 9; Day-1 scope: just look for the file)
    holder_path = Path(f"{prefix}.lock.holder")
    if holder_path.exists():
        try:
            holder_pid = holder_path.read_text().strip() or "unknown"
        except OSError:
            holder_pid = "unknown"
        return CheckResult(
            "output prefix lock free",
            "WARN",
            f"another process (PID {holder_pid}) currently holds the lock; "
            "extract would block until released",
        )

    candidates = [
        Path(f"{prefix}.geno"),
        Path(f"{prefix}.snp"),
        Path(f"{prefix}.ind"),
        Path(f"{prefix}.pseudohaploid.json"),
    ]
    existing = [p for p in candidates if p.exists()]
    if existing:
        return CheckResult(
            "output prefix free",
            "WARN",
            f"{len(existing)} output file(s) exist at prefix; extract would refuse without --overwrite",
        )

    return CheckResult("output prefix writable + free", "PASS", "")


def _emit_results(results: list[CheckResult], *, json_output: bool) -> int:
    any_fail = any(r.status == "FAIL" for r in results)
    if json_output:
        payload = {
            "checks": [
                {"name": r.name, "status": r.status, "detail": r.detail}
                for r in results
            ],
            "exit_code": 1 if any_fail else 0,
        }
        json.dump(payload, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        sys.stdout.write("status\tcheck\tdetail\n")
        for r in results:
            sys.stdout.write(f"{r.status}\t{r.name}\t{r.detail}\n")
    return 1 if any_fail else 0


__all__ = ["CheckResult", "run_validate"]
