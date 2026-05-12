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

import hashlib
import importlib.resources as res
import logging
import os
import re
import time
from pathlib import Path
from typing import Final

import pysam

from .counters import Stage1InputFilters, Stage1LiftCounters
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
        input_filters={
            "palindrome_drops": input_filter_counters.palindrome_drops,
            "non_snp_drops": input_filter_counters.non_snp_drops,
            "non_autosome_drops": input_filter_counters.non_autosome_drops,
        },
    )


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
    "chain_file_path",
    "get_bundled_chain_path",
    "lift_aadr_sites",
    "parse_picard_stderr",
    "parse_rejected_vcf",
    "resolve_chain_for_extract",
]
