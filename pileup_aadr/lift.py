"""Stage 1: Picard LiftoverVcf + chain-file resolution.

Day 2 scope: chain-file resolution helpers + bundled-chain SHA verification.
The full `lift_aadr_sites` orchestration lands on Day 3 (after sites_vcf.py).
"""
from __future__ import annotations

import hashlib
import importlib.resources as res
import logging
import os
from pathlib import Path

from .errors import ChainFileNotFound, ChainFileSHAError

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


__all__ = [
    "chain_file_path",
    "get_bundled_chain_path",
    "resolve_chain_for_extract",
]
