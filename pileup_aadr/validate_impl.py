"""`validate` subcommand implementation — 10 pre-flight checks per LLD §16.

Day-2 scope: tool-version probes via `tool_wrapper` + ref-FASTA findability + build
verification via `ref_resolve` + chain-file resolution via `lift`. The stubbed
binary-presence checks from Day 1 have been replaced with real version checks.
"""
from __future__ import annotations

import json
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import pandas as pd

from .errors import (
    PileupAadrError,
    ToolNotFoundError,
    ToolVersionError,
)
from .format_detect import (
    BuildOverride,
    detect_aadr_build,
    detect_bam_build,
    detect_bam_format,
    parse_aadr_snp,
)
from .lift import resolve_chain_for_extract
from .ref_resolve import verify_fasta_matches_bam_build
from .tool_wrapper import (
    JAVA_SPEC,
    PICARD_SPEC,
    PILEUPCALLER_SPEC,
    SAMTOOLS_SPEC,
    ToolSpec,
    ToolWrapper,
)

log = logging.getLogger(__name__)


@dataclass
class CheckResult:
    """One pre-flight check result.

    `resolved` carries an artifact path back to follow-up checks (currently
    only the resolved target FASTA Path, used by `_check_target_fasta_dict`).
    Excluded from emit_results output so the user-facing TSV/JSON shapes
    don't change.
    """

    name: str
    status: Literal["PASS", "WARN", "FAIL", "SKIP"]
    detail: str = ""
    resolved: Path | None = None


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

    # --- Tool dependencies (Day-2 fix: real version probes via tool_wrapper) ---
    results.append(_check_tool(SAMTOOLS_SPEC))
    if no_lift:
        results.append(CheckResult("picard binary", "SKIP", "no-lift fast path"))
        results.append(CheckResult("java binary", "SKIP", "no-lift fast path"))
    else:
        results.append(_check_tool(PICARD_SPEC))
        results.append(_check_tool(JAVA_SPEC))
    results.append(_check_tool(PILEUPCALLER_SPEC))

    # --- Tool flag probes (v0.2 #1: catch tool-CLI drift at validate time) ---
    # Run each tool's --help and confirm the specific flags pileup-aadr will
    # use are present. Catches the bug class where samtools/Picard/etc. rename
    # or drop a flag in a new release (e.g., issue #2: `samtools mpileup -@`
    # was always invalid; we shipped it for two patch releases anyway).
    results.append(_check_tool_flags_samtools_mpileup())
    results.append(_check_tool_flags_pileupcaller())
    results.append(_check_tool_flags_mosdepth())
    if no_lift:
        results.append(
            CheckResult("picard LiftoverVcf flag probe", "SKIP", "no-lift fast path"),
        )
    else:
        results.append(_check_tool_flags_picard_liftover())

    # --- Chain (Day-2: bundled-chain SHA verification + 3-tier resolution) ---
    results.append(_check_chain(chain_path))

    # --- Target FASTA (Day-2: real build verification via ref_resolve) ---
    fasta_check = _check_ref_fasta(ref_fasta, bam, bam_build)
    results.append(fasta_check)

    # --- Target FASTA .dict (skipped on no-lift fast path; Picard isn't called) ---
    results.append(_check_target_fasta_dict(fasta_check, no_lift))

    # --- Output prefix ---
    if output_prefix is None:
        results.append(CheckResult("output prefix", "SKIP", "no -o passed"))
    else:
        results.append(_check_output_prefix(output_prefix))

    return _emit_results(results, json_output=json_output)


def _aadr_parse_pass(aadr_df: object) -> CheckResult:
    """parse_aadr_snp succeeded — render a PASS with summary."""
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


def _probe_help_for_flags(
    name: str,
    invocation: list[str],
    required_flags: tuple[str, ...],
    *,
    timeout: float = 10.0,
) -> CheckResult:
    """Run `invocation` (the help command), grep stdout+stderr for each
    required flag. Returns PASS if all flags found, FAIL listing the missing
    ones. Tool exit code is ignored — many tools exit non-zero on --help.

    `name` is the human-readable check name shown in the report.
    `required_flags` are matched as whole-word substrings (e.g., '-@'
    matches '-@' or '-@,' but not 'foo-@bar').
    """
    import re
    import subprocess
    try:
        proc = subprocess.run(
            invocation, capture_output=True, text=True,
            timeout=timeout, check=False,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return CheckResult(name, "FAIL", f"could not run {invocation[0]}: {e}")
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    missing: list[str] = []
    for flag in required_flags:
        # Word-boundary-ish match: flag preceded/followed by non-alnum or line edge.
        pattern = rf"(?:^|[^A-Za-z0-9_-]){re.escape(flag)}(?:$|[^A-Za-z0-9_-])"
        if not re.search(pattern, combined):
            missing.append(flag)
    if missing:
        return CheckResult(
            name, "FAIL",
            f"flags not found in --help output: {', '.join(missing)}",
        )
    return CheckResult(
        name, "PASS",
        f"all {len(required_flags)} flag(s) present: {', '.join(required_flags)}",
    )


def _check_tool_flags_samtools_mpileup() -> CheckResult:
    """Verify samtools mpileup accepts -B, -q, -Q, -R, -f, -l (the flags
    pileup-aadr's Stage 3 actually uses). Issue #2 motivation: `mpileup -@`
    was rejected silently for two releases."""
    return _probe_help_for_flags(
        name="samtools mpileup flag probe",
        invocation=["samtools", "mpileup", "--help"],
        required_flags=("-B", "-q", "-Q", "-R", "-f", "-l"),
    )


def _check_tool_flags_pileupcaller() -> CheckResult:
    """Verify pileupCaller accepts the flags Stage 3 uses."""
    return _probe_help_for_flags(
        name="pileupCaller flag probe",
        invocation=["pileupCaller", "--help"],
        required_flags=(
            "--randomDiploid", "--seed", "--sampleNames",
            "--samplePopName", "-e",
        ),
    )


def _check_tool_flags_mosdepth() -> CheckResult:
    """Verify mosdepth accepts the flags `coverage` subcommand uses."""
    return _probe_help_for_flags(
        name="mosdepth flag probe",
        invocation=["mosdepth", "--help"],
        required_flags=("--threads", "--no-per-base", "--quantize", "--by"),
    )


def _check_tool_flags_picard_liftover() -> CheckResult:
    """Verify Picard LiftoverVcf accepts the flags Stage 1 uses.

    Picard subcommand flag probing requires the JAR + java; if either is
    unresolvable, the check FAILs with a pointer at the missing dependency.
    """
    import shutil

    from .tool_wrapper import _resolve_picard_jar

    if shutil.which("java") is None:
        return CheckResult(
            "picard LiftoverVcf flag probe", "FAIL",
            "java not on PATH (needed to invoke picard.jar)",
        )
    try:
        jar = _resolve_picard_jar()
    except PileupAadrError as e:
        return CheckResult("picard LiftoverVcf flag probe", "FAIL", e.why)

    return _probe_help_for_flags(
        name="picard LiftoverVcf flag probe",
        invocation=["java", "-jar", str(jar), "LiftoverVcf", "--help"],
        required_flags=(
            "--INPUT", "--OUTPUT", "--CHAIN", "--REFERENCE_SEQUENCE",
            "--REJECT", "--RECOVER_SWAPPED_REF_ALT",
            "--MAX_RECORDS_IN_RAM",
        ),
        timeout=30.0,  # Picard JVM startup is slow
    )


def _check_tool(spec: ToolSpec) -> CheckResult:
    """Probe a tool's binary + version via ToolWrapper (Day-2 fix: real version check)."""
    try:
        observed = ToolWrapper(spec).version()
        return CheckResult(
            f"{spec.binary} version",
            "PASS",
            f"{observed} (>= {spec.min_version} required; tested against {spec.tested_against})",
        )
    except (ToolNotFoundError, ToolVersionError, PileupAadrError) as e:
        return CheckResult(f"{spec.binary} version", "FAIL", e.why)


def _check_chain(cli_chain: Path | None) -> CheckResult:
    """Verify the chain file resolves cleanly (bundled or user-supplied) + SHA OK."""
    try:
        resolved = resolve_chain_for_extract(cli_chain=cli_chain)
        if cli_chain is not None:
            return CheckResult(
                "chain file", "PASS", f"user-supplied: {resolved} (SHA not enforced)"
            )
        return CheckResult(
            "chain file", "PASS", f"bundled: {resolved.name} (SHA verified)"
        )
    except PileupAadrError as e:
        return CheckResult("chain file", "FAIL", e.why)


def _check_ref_fasta(
    cli_ref: Path | None,
    bam: Path | None,
    bam_build: str | None,
) -> CheckResult:
    """Resolve target FASTA + verify chr1 length matches BAM build (Day-2 fix).

    For pre-flight purposes we accept any of the three resolution paths (--ref-fasta,
    env, BAM @PG) and report which one resolved. If `bam_build` is unknown (BAM detection
    failed), we skip the build-match verification but still report whether a file exists.

    The resolved Path is also attached to `CheckResult.resolved` for the
    follow-up `_check_target_fasta_dict` call.
    """
    if cli_ref is not None:
        if not cli_ref.exists():
            return CheckResult(
                "ref FASTA findable", "FAIL", f"--ref-fasta {cli_ref} not found"
            )
        if bam_build is not None:
            try:
                verify_fasta_matches_bam_build(cli_ref, bam_build)  # type: ignore[arg-type]
                return CheckResult(
                    "ref FASTA findable + build match", "PASS", str(cli_ref),
                    resolved=cli_ref,
                )
            except PileupAadrError as e:
                return CheckResult("ref FASTA build mismatch", "FAIL", e.why)
        return CheckResult(
            "ref FASTA findable",
            "WARN",
            f"--ref-fasta: {cli_ref} (build verification skipped — BAM build unknown)",
            resolved=cli_ref,
        )

    env_dir = os.environ.get("PILEUP_AADR_REF_DIR")
    if env_dir and bam_build is not None:
        candidate = Path(env_dir) / f"{bam_build}.fa"
        if candidate.exists():
            try:
                verify_fasta_matches_bam_build(candidate, bam_build)  # type: ignore[arg-type]
                return CheckResult(
                    "ref FASTA findable + build match", "PASS",
                    f"$PILEUP_AADR_REF_DIR: {candidate}",
                    resolved=candidate,
                )
            except PileupAadrError as e:
                return CheckResult("ref FASTA build mismatch", "FAIL", e.why)

    # BAM @PG fallback path — only attempt if we have a BAM to inspect.
    # The full extraction logic lives in ref_resolve._extract_ref_from_bam_pg
    # but for validate's pre-flight we do a lightweight inspection that
    # doesn't raise on missing path.
    if bam is not None and bam_build is not None:
        from .ref_resolve import _extract_ref_from_bam_pg

        try:
            extracted = _extract_ref_from_bam_pg(bam)
            try:
                verify_fasta_matches_bam_build(extracted, bam_build)  # type: ignore[arg-type]
                return CheckResult(
                    "ref FASTA findable + build match", "PASS",
                    f"BAM @PG: {extracted}",
                    resolved=extracted,
                )
            except PileupAadrError as e:
                return CheckResult("ref FASTA build mismatch", "FAIL", e.why)
        except PileupAadrError as e:
            return CheckResult("ref FASTA findable", "WARN", e.why)

    return CheckResult(
        "ref FASTA findable",
        "WARN",
        "no --ref-fasta, no $PILEUP_AADR_REF_DIR, and BAM @PG yielded no candidates",
    )


def _check_target_fasta_dict(
    fasta_check: CheckResult,
    no_lift: bool,
) -> CheckResult:
    """Verify (don't auto-generate) the .dict alongside the resolved FASTA.

    On the no-lift fast path, Picard isn't invoked so .dict isn't needed.
    Otherwise: PASS if a .dict exists; WARN if missing-but-FASTA-dir-writable
    (extract will auto-generate it); FAIL if the FASTA itself failed to
    resolve in the prior check.
    """
    if no_lift:
        return CheckResult(
            "target FASTA .dict",
            "SKIP",
            "no-lift fast path (Picard not invoked; .dict not needed)",
        )
    if fasta_check.resolved is None:
        return CheckResult(
            "target FASTA .dict",
            "FAIL",
            "FASTA did not resolve in prior check; cannot check .dict",
        )

    from .dict_resolve import find_existing_dict, find_or_user_cache_dict_path

    fasta = fasta_check.resolved
    existing = find_existing_dict(fasta)
    if existing is not None:
        return CheckResult("target FASTA .dict", "PASS", str(existing))

    # No existing .dict — would extract auto-generate? Yes if a write target exists.
    target = find_or_user_cache_dict_path(fasta)
    return CheckResult(
        "target FASTA .dict",
        "WARN",
        (
            f"no .dict alongside {fasta}; extract will auto-generate at {target} "
            "(one-time ~23s on hg38)"
        ),
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
