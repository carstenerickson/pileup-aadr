"""Tests for lift.py — chain-file resolution + bundled-chain SHA verification.

Maps to LLD test #39 (bundled chain SHA matches sidecar) + #42 (chain-resolution
flag matrix).
"""
from __future__ import annotations

import hashlib
import importlib.resources as res
from pathlib import Path

import pytest

from pileup_aadr.errors import ChainFileNotFound, ChainFileSHAError
from pileup_aadr.lift import (
    chain_file_path,
    get_bundled_chain_path,
    resolve_chain_for_extract,
)

# --- Bundled chain integrity (LLD test #39) ---


def test_bundled_chain_sha_matches_sidecar() -> None:
    """SHA in the .sha256 sidecar matches the chain file's actual SHA."""
    chain_resource = res.files("pileup_aadr") / "data" / "hg19ToHg38.over.chain.gz"
    sha_resource = res.files("pileup_aadr") / "data" / "hg19ToHg38.over.chain.gz.sha256"
    expected_sha = sha_resource.read_text().strip().split()[0]
    actual_sha = hashlib.sha256(chain_resource.read_bytes()).hexdigest()
    assert actual_sha == expected_sha, (
        f"Bundled chain SHA mismatch:\n  expected: {expected_sha}\n  actual:   {actual_sha}\n"
        "Run `make refresh-chain-sha` to regenerate the sidecar."
    )


def test_bundled_chain_canonical_size() -> None:
    """LLD: UCSC's hg19ToHg38.over.chain.gz has been ~223 KB for a decade.

    Per outside-review F4 fix: use len(read_bytes()) instead of .stat().st_size since
    Traversable.stat() is not part of importlib.resources contract (breaks for
    zip-archived packages).
    """
    chain_bytes = (res.files("pileup_aadr") / "data" / "hg19ToHg38.over.chain.gz").read_bytes()
    assert 200_000 <= len(chain_bytes) <= 250_000, (
        f"Chain file size {len(chain_bytes)} outside expected 200-250 KB range; "
        "verify this is the correct UCSC chain"
    )


def test_get_bundled_chain_path_succeeds() -> None:
    """Default startup verification passes on the committed bundled chain."""
    p = get_bundled_chain_path()
    assert p.exists()
    assert p.name == "hg19ToHg38.over.chain.gz"


# --- 3-tier resolution (LLD test #42 — chain-resolution flag matrix) ---


def test_chain_file_path_returns_bundled_when_no_overrides() -> None:
    """No --chain, no $PILEUP_AADR_CHAIN_DIR → bundled chain returned."""
    result = chain_file_path(cli_chain=None, env_chain_dir=None)
    assert result.name == "hg19ToHg38.over.chain.gz"


def test_chain_file_path_uses_cli_chain_when_given(tmp_path: Path) -> None:
    """--chain PATH explicit override beats env + bundled."""
    user_chain = tmp_path / "custom.chain.gz"
    user_chain.write_bytes(b"fake chain bytes for testing")
    result = chain_file_path(cli_chain=user_chain, env_chain_dir=None)
    assert result == user_chain


def test_chain_file_path_uses_env_when_no_cli(tmp_path: Path) -> None:
    """$PILEUP_AADR_CHAIN_DIR resolved when --chain not set."""
    env_dir = tmp_path / "chains"
    env_dir.mkdir()
    env_chain = env_dir / "hg19ToHg38.over.chain.gz"
    env_chain.write_bytes(b"fake env chain")
    result = chain_file_path(cli_chain=None, env_chain_dir=env_dir)
    assert result == env_chain


def test_chain_file_path_falls_through_to_bundled_when_env_missing_file(
    tmp_path: Path,
) -> None:
    """$PILEUP_AADR_CHAIN_DIR set but the file isn't there → fall through to bundled."""
    env_dir = tmp_path / "chains_missing"
    env_dir.mkdir()
    # Don't create the chain file in env_dir
    result = chain_file_path(cli_chain=None, env_chain_dir=env_dir)
    assert result.name == "hg19ToHg38.over.chain.gz"
    # Should be the bundled one (from the package), not the env one
    assert "pileup_aadr/data" in str(result)


def test_chain_file_path_raises_chain_not_found_for_bad_cli_path(tmp_path: Path) -> None:
    """--chain pointing at a nonexistent path raises ChainFileNotFound."""
    bogus = tmp_path / "does_not_exist.chain.gz"
    with pytest.raises(ChainFileNotFound, match="does not exist"):
        chain_file_path(cli_chain=bogus, env_chain_dir=None)


def test_chain_file_path_strict_sha_with_user_chain_mismatch(tmp_path: Path) -> None:
    """--strict-chain-sha + user-supplied chain that doesn't match → ChainFileSHAError."""
    user_chain = tmp_path / "custom.chain.gz"
    user_chain.write_bytes(b"definitely not the canonical UCSC chain bytes")
    with pytest.raises(ChainFileSHAError, match="strict-chain-sha"):
        chain_file_path(cli_chain=user_chain, env_chain_dir=None, strict_sha=True)


def test_chain_file_path_insecure_skips_strict_check(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """--insecure-chain takes precedence over --strict-chain-sha (skips, with warning)."""
    user_chain = tmp_path / "custom.chain.gz"
    user_chain.write_bytes(b"not the canonical chain")
    import logging
    caplog.set_level(logging.WARNING, logger="pileup_aadr.lift")
    result = chain_file_path(cli_chain=user_chain, strict_sha=True, insecure=True, env_chain_dir=None)
    assert result == user_chain
    assert any("verification skipped" in r.message for r in caplog.records)


# --- env-var resolution helper ---


def test_resolve_chain_for_extract_reads_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """resolve_chain_for_extract pulls $PILEUP_AADR_CHAIN_DIR from os.environ."""
    env_dir = tmp_path / "chains"
    env_dir.mkdir()
    env_chain = env_dir / "hg19ToHg38.over.chain.gz"
    env_chain.write_bytes(b"env chain bytes")
    monkeypatch.setenv("PILEUP_AADR_CHAIN_DIR", str(env_dir))
    result = resolve_chain_for_extract(cli_chain=None)
    assert result == env_chain


def test_resolve_chain_for_extract_no_env_no_cli_returns_bundled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both env and CLI unset → bundled chain returned."""
    monkeypatch.delenv("PILEUP_AADR_CHAIN_DIR", raising=False)
    result = resolve_chain_for_extract(cli_chain=None)
    assert result.name == "hg19ToHg38.over.chain.gz"
