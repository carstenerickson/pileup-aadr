"""Fixtures for the real-binary integration tests.

Skipped automatically (and silently) if the required external binary is missing
from PATH or PICARD_JAR is unset. CI's bio-tools job installs all four via
conda + sets PICARD_JAR; local runs `pip install -e '.[dev]'` skip these.

The chr22 FASTA fixture downloads UCSC's hg38 chr22 (~14 MB compressed) on
first use and caches it under `tests/integration/fixtures/`. The download is
~1s on CI runners; subsequent runs are a no-op.
"""
from __future__ import annotations

import shutil
import subprocess
import urllib.request
from pathlib import Path

import pytest

INTEGRATION_DIR = Path(__file__).parent
INTEGRATION_FIXTURES = INTEGRATION_DIR / "fixtures"

_HG38_CHR22_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"


def _has_binary(name: str) -> bool:
    return shutil.which(name) is not None


def _has_picard() -> bool:
    """Picard is detected by attempting the same resolver `tool_wrapper` uses.

    A `picard` shell wrapper on PATH is necessary but not sufficient — Homebrew's
    `picard` wrapper exists without a colocated JAR, which would skip the test
    for the wrong reason. Run the actual resolver and return True only if it
    yields a real JAR path.
    """
    from pileup_aadr.tool_wrapper import _resolve_picard_jar
    try:
        _resolve_picard_jar()
    except Exception:
        return False
    return True


requires_samtools = pytest.mark.skipif(
    not _has_binary("samtools"), reason="samtools not on PATH",
)
requires_pileupcaller = pytest.mark.skipif(
    not _has_binary("pileupCaller"), reason="pileupCaller not on PATH",
)
requires_mosdepth = pytest.mark.skipif(
    not _has_binary("mosdepth"), reason="mosdepth not on PATH",
)
requires_picard = pytest.mark.skipif(
    not _has_picard(),
    reason="Picard not findable (set PICARD_JAR or install via bioconda)",
)


@pytest.fixture(scope="session")
def hg38_chr22_fasta() -> Path:
    """Download UCSC's hg38 chr22 FASTA + faidx + create .dict.

    Cached under `tests/integration/fixtures/chr22.fa` after the first run.
    The download is ~14 MB compressed, ~50 MB uncompressed; faidx +
    CreateSequenceDictionary together take ~1s on chr22.
    """
    INTEGRATION_FIXTURES.mkdir(parents=True, exist_ok=True)
    fasta = INTEGRATION_FIXTURES / "chr22.fa"
    if not fasta.exists():
        gz = INTEGRATION_FIXTURES / "chr22.fa.gz"
        urllib.request.urlretrieve(_HG38_CHR22_URL, gz)
        subprocess.run(["gunzip", "-f", str(gz)], check=True)
    fai = Path(f"{fasta}.fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    return fasta


@pytest.fixture(scope="session")
def hg38_chr22_bam(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A minimal hg38-headered BAM with chr1+chr22 @SQ entries (no reads).

    Sufficient for pileup-aadr's pre-flight + Stage 1/2/4 pipeline; Stage 3
    pileupCaller will return 0 non-missing calls (no reads in the BAM) but
    still emits a parseable summary stderr block.
    """
    import pysam

    out = tmp_path_factory.mktemp("integration") / "user.bam"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [
            {"SN": "chr1", "LN": 248_956_422},
            {"SN": "chr22", "LN": 50_818_468},  # hg38 chr22 length
        ],
        "RG": [{"ID": "1", "SM": "TestSample", "LB": "lib1", "PL": "ILLUMINA"}],
    }
    with pysam.AlignmentFile(str(out), "wb", header=header):
        pass
    pysam.index(str(out))
    return out


@pytest.fixture(scope="session")
def aadr_chr22_slice() -> Path:
    """Reuse the unit-test fixture: 50 chr22 sites in hg19 coords."""
    return Path(__file__).parent.parent / "fixtures" / "aadr_chr22_slice.snp"
