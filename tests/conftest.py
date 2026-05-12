"""Session fixtures for pytest. Most live under tests/fixtures/.

Day 1 scope: only the AADR slice fixture is created, since the only modules
implemented today are format_detect / inspect / validate (file-side checks).
The synthetic BAM + chain fixtures are spec'd in LLD §18 B5/B6 and land on
later days when the corresponding stage modules are implemented.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture(autouse=True)
def quiet_logging() -> None:
    """Silence INFO logging during tests. Without this, format_detect.parse_aadr_snp
    emits 'Loaded AADR .snp: N rows...' to stderr, which CliRunner mixes into stdout
    and breaks JSON parsing in the smoke tests. Tests that explicitly want INFO can
    re-configure via configure_logging(level=logging.INFO) inside the test body.
    """
    from pileup_aadr.logging_config import configure_logging

    configure_logging(level=logging.WARNING)


@pytest.fixture(scope="session")
def aadr_chr22_slice() -> Path:
    """A small synthetic AADR .snp slice for Day-1 inspect/validate tests.

    Format matches AADR v66 1240k:
        whitespace-padded 6 columns: rsid, chrom_int, gen_morgans, pos_bp, REF, ALT
        chrom encoding: 1-22 numeric, 23=X, 24=Y, 90=mt
    """
    return FIXTURES_DIR / "aadr_chr22_slice.snp"


@pytest.fixture
def tmp_aadr_snp(tmp_path: Path, aadr_chr22_slice: Path) -> Path:
    """Copy the slice into a writable tmp_path so tests can mutate without touching the fixture."""
    out = tmp_path / "aadr_test.snp"
    out.write_bytes(aadr_chr22_slice.read_bytes())
    return out
