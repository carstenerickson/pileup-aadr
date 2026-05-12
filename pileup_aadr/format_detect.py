"""Pre-stage detection of input formats + builds.

Cheapest, fastest module — pure metadata reads, no genotype I/O. Called once at
orchestrator startup before any expensive work begins. Failures surface immediately
with clear errors per HLD §"Error class taxonomy".
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Final, Literal

import pandas as pd
import pysam

from .errors import (
    AADRDuplicateRsidError,
    AADRParseError,
    BAMParseError,
    BAMSampleNameAmbiguous,
    UnsupportedAADRBuild,
    UnsupportedReferenceBuild,
)

log = logging.getLogger(__name__)

BuildOverride = Literal["hg19", "hg38", "auto"]

# Canonical chr1 lengths per assembly (hardcoded; chromosome lengths are byte-stable)
HG19_CHR1_LENGTH: Final[int] = 249_250_621
HG38_CHR1_LENGTH: Final[int] = 248_956_422
_BUILD_TOLERANCE_BP: Final[int] = 1_000_000  # ±1Mb tolerance for patch-level differences

# IID character set restriction for EIGENSOFT compatibility
_SAFE_IID_RE: Final[re.Pattern[str]] = re.compile(r"[^A-Za-z0-9_.-]")

# Chromosome normalization: maps any input form to canonical chrN.
# Returns None for unrecognized inputs (alt-contigs, decoys, etc.).
_CHROM_NORMALIZE: Final[dict[str, str]] = {
    **{str(i): f"chr{i}" for i in range(1, 23)},
    **{f"chr{i}": f"chr{i}" for i in range(1, 23)},
    "23": "chrX", "X": "chrX", "chrX": "chrX",
    "24": "chrY", "Y": "chrY", "chrY": "chrY",
    "90": "chrM", "MT": "chrM", "M": "chrM", "chrM": "chrM", "chrMT": "chrM",
    "91": "chrXY",  # PAR; rare
}


def normalize_chrom(raw: str) -> str | None:
    """Normalize a chromosome name to its canonical chrN form.

    Args:
        raw: chromosome name in any common encoding ("1", "chr1", "23", "X", "MT", etc.)

    Returns:
        Canonical chrN form (e.g., "chr1", "chrX", "chrM"), or None for unrecognized
        inputs (alt-contigs, decoys). Returning None lets callers filter alt-contigs
        cleanly instead of carrying them through.
    """
    return _CHROM_NORMALIZE.get(raw)


def detect_bam_format(bam_path: Path) -> Literal["BAM", "CRAM"]:
    """Detect BAM vs CRAM from magic bytes. File extension is a hint, not authoritative.

    Args:
        bam_path: path to a candidate alignment file

    Returns:
        "BAM" or "CRAM" string

    Raises:
        BAMParseError: file is neither BAM nor CRAM (e.g., truncated, wrong format)
    """
    with open(bam_path, "rb") as f:
        magic = f.read(4)
    if magic == b"\x1f\x8b\x08\x04":
        # gzip magic with extra-field flag — could be BAM (always has extra field) or vanilla gzip.
        # BAM-specific: opens cleanly under "rb" mode and `is_bam` returns True.
        with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
            if bam.is_bam:
                return "BAM"
            raise BAMParseError(
                what=str(bam_path),
                why="gzip-compressed but not a valid BAM",
                fix="Verify file integrity with `samtools quickcheck`",
            )
    if magic == b"CRAM":
        return "CRAM"
    raise BAMParseError(
        what=str(bam_path),
        why=f"unrecognized magic bytes {magic!r}",
        fix="Expected BAM (gzip-compressed) or CRAM; verify file is an aligned-reads file",
    )


def detect_bam_build(
    bam_path: Path, override: BuildOverride = "auto"
) -> Literal["hg19", "hg38"]:
    """Detect BAM coordinate system from @SQ chr1 length.

    Args:
        bam_path: aligned BAM/CRAM with @SQ headers
        override: "auto" (detect), "hg19", or "hg38" (skip detection)

    Returns:
        "hg19" or "hg38"

    Raises:
        UnsupportedReferenceBuild: chr1 length matches neither hg19 nor hg38 ± 1Mb
            (e.g., T2T-CHM13). Pass --bam-build explicitly to override detection if
            the BAM is hg19/hg38-compatible despite differing chr1 length.
    """
    if override != "auto":
        return override  # type: ignore[return-value]

    # v2.2 M13 fix: "r" auto-detects BAM vs CRAM; "rb" was BAM-only and crashed on CRAM.
    # Header reads don't need reference_filename even for CRAM (header is in the file).
    with pysam.AlignmentFile(str(bam_path), "r", check_sq=False) as bam:
        sq_records = bam.header.get("SQ", [])
        chr1_record = next(
            (sq for sq in sq_records if sq["SN"] in ("chr1", "1")),
            None,
        )
        if chr1_record is None:
            raise UnsupportedReferenceBuild(
                what=str(bam_path),
                why="BAM @SQ has no chr1 / 1 record; cannot determine build",
                fix="Pass --bam-build hg19|hg38 to override detection",
            )
        chr1_length = int(chr1_record["LN"])

    if abs(chr1_length - HG19_CHR1_LENGTH) <= _BUILD_TOLERANCE_BP:
        return "hg19"
    if abs(chr1_length - HG38_CHR1_LENGTH) <= _BUILD_TOLERANCE_BP:
        return "hg38"
    raise UnsupportedReferenceBuild(
        what=str(bam_path),
        why=(
            f"BAM @SQ chr1 length {chr1_length} matches neither hg19 ({HG19_CHR1_LENGTH}) "
            f"nor hg38 ({HG38_CHR1_LENGTH}) within ±1 Mb"
        ),
        fix=(
            "Pass --bam-build hg19|hg38 to override (if the BAM is hg19/hg38-compatible "
            "despite differing chr1 length, e.g., T2T-CHM13 fragments)"
        ),
    )


def detect_bam_sample_name(bam_path: Path, explicit: str | None = None) -> str:
    """Resolve sample name from --sample-name → @RG SM: → filename stem.

    Args:
        bam_path: aligned BAM/CRAM
        explicit: --sample-name CLI value (None if not given)

    Returns:
        Sample name suitable for EIGENSTRAT .ind IID column. Characters outside
        [A-Za-z0-9_.-] are replaced with '_' and a stderr WARNING is logged.

    Raises:
        BAMSampleNameAmbiguous: multiple @RG SM: values disagree AND --sample-name not given
    """
    if explicit is not None:
        # Explicit override always wins; emit INFO line if it differs from @RG SM:
        rg_sms = _extract_rg_sms(bam_path)
        if rg_sms and explicit not in rg_sms:
            log.info(
                "Using --sample-name '%s' (overriding BAM @RG SM: %s)",
                explicit, sorted(rg_sms),
            )
        return _sanitize_iid(explicit)

    rg_sms = _extract_rg_sms(bam_path)
    if len(rg_sms) == 1:
        return _sanitize_iid(next(iter(rg_sms)))
    if len(rg_sms) > 1:
        raise BAMSampleNameAmbiguous(
            what=str(bam_path),
            why=f"BAM has {len(rg_sms)} disagreeing @RG SM: values: {sorted(rg_sms)}",
            fix="Pass --sample-name explicitly to disambiguate",
        )
    # No @RG SM at all → fall back to filename stem
    return _sanitize_iid(bam_path.stem)


def _extract_rg_sms(bam_path: Path) -> set[str]:
    """Return the set of distinct SM values across all @RG lines (empty set if none).

    Uses "r" auto-detect mode (M13 fix) to support BAM + CRAM uniformly.
    """
    with pysam.AlignmentFile(str(bam_path), "r", check_sq=False) as bam:
        rgs = bam.header.get("RG", [])
    return {rg["SM"] for rg in rgs if "SM" in rg}


def _sanitize_iid(raw: str) -> str:
    """Replace IID characters outside [A-Za-z0-9_.-] with '_' and warn on substitutions."""
    cleaned = _SAFE_IID_RE.sub("_", raw)
    if cleaned != raw:
        log.warning(
            "IID '%s' contains characters outside [A-Za-z0-9_.-]; sanitized to '%s'",
            raw, cleaned,
        )
    return cleaned


_AADR_COLUMNS: Final[list[str]] = ["rsid", "chrom_int", "gen_morgans", "pos_bp", "ref", "alt"]


def parse_aadr_snp(aadr_snp_path: Path) -> pd.DataFrame:
    """Load AADR .snp into a DataFrame indexed by rsid.

    Format (verified empirically v2.1):
        6 columns whitespace-separated (multi-space-padded for column alignment, NOT
        tab-separated). No header. Lines starting with '#' tolerated as comments.

    Args:
        aadr_snp_path: AADR .snp file

    Returns:
        DataFrame with 6 typed columns + rsid index. dtypes:
            rsid: str (index)
            chrom_int: str (kept as string to handle "1"-"22","23","24","90","91")
            gen_morgans: float
            pos_bp: int
            ref: str (single ACGT char)
            alt: str (single ACGT char)

    Raises:
        AADRParseError: row has wrong column count, non-ACGT alleles, or unparseable position
        AADRDuplicateRsidError: duplicate rsid in file (Stage 4 join requires unique IDs)
    """
    rows: list[tuple[str, str, float, int, str, str]] = []
    with open(aadr_snp_path) as f:
        for lineno, line in enumerate(f, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()  # str.split() handles multi-space + tabs
            if len(parts) != 6:
                raise AADRParseError(
                    what=f"{aadr_snp_path}:{lineno}",
                    why=f"expected 6 columns, got {len(parts)}: {stripped[:80]!r}",
                    fix="Verify AADR file integrity; re-download from Reich Lab if needed",
                )
            rsid, chrom_int, gen_str, pos_str, ref, alt = parts
            try:
                gen = float(gen_str)
                pos = int(pos_str)
            except ValueError as e:
                raise AADRParseError(
                    what=f"{aadr_snp_path}:{lineno}",
                    why=f"unparseable genetic-distance or position: {e}",
                    fix=(
                        "Verify AADR file integrity; columns 3 (Morgans) and 4 (bp) "
                        "must be numeric"
                    ),
                ) from e
            if ref not in "ACGT" or alt not in "ACGT":
                raise AADRParseError(
                    what=f"{aadr_snp_path}:{lineno}",
                    why=f"non-ACGT alleles: REF={ref!r} ALT={alt!r}",
                    fix="AADR .snp must contain only biallelic SNPs with ACGT alleles",
                )
            rows.append((rsid, chrom_int, gen, pos, ref, alt))

    df = pd.DataFrame(rows, columns=_AADR_COLUMNS)

    # Duplicate-rsID startup invariant (matches pgen-samplebind's input[0] uniqueness check)
    dup_mask = df["rsid"].duplicated(keep=False)
    if dup_mask.any():
        first_dup = df.loc[dup_mask, "rsid"].iloc[0]
        first_indices = df.index[df["rsid"] == first_dup].tolist()
        raise AADRDuplicateRsidError(
            what=f"{aadr_snp_path} (rsid={first_dup})",
            why=(
                f"rsid '{first_dup}' appears at rows "
                f"{first_indices[0] + 1} and {first_indices[1] + 1}"
            ),
            fix="Verify AADR file integrity; AADR .snp must have unique SNP names",
        )
    df = df.set_index("rsid", verify_integrity=True)
    log.info(
        "Loaded AADR .snp: %d rows, %d unique chromosomes",
        len(df), df["chrom_int"].nunique(),
    )
    return df


def detect_aadr_build(
    aadr_df: pd.DataFrame, override: BuildOverride = "auto"
) -> Literal["hg19", "hg38"]:
    """Detect AADR coordinate system from chr1 max bp position.

    Args:
        aadr_df: DataFrame from `parse_aadr_snp`
        override: "auto" (detect from data), "hg19", or "hg38" (skip detection)

    Returns:
        "hg19" or "hg38"

    Raises:
        UnsupportedAADRBuild: chr1 max position matches neither hg19 nor hg38 within tolerance,
            or there are no chr1 rows.
    """
    if override != "auto":
        return override  # type: ignore[return-value]

    chr1_rows = aadr_df[aadr_df["chrom_int"] == "1"]
    if chr1_rows.empty:
        raise UnsupportedAADRBuild(
            what="(AADR DataFrame)",
            why="no chr1 (chrom_int='1') rows found; cannot determine build",
            fix="Pass --aadr-build hg19|hg38 to override detection",
        )
    chr1_max = int(chr1_rows["pos_bp"].max())
    # The hg19 and hg38 chr1 windows overlap heavily (hg19 = 249,250,621; hg38 = 248,956,422;
    # difference 294 KB). A simple "<= length AND >= length - tolerance" check picks hg19
    # for any position in the overlap. Pick the closest assembly instead.
    dist_hg19 = abs(chr1_max - HG19_CHR1_LENGTH)
    dist_hg38 = abs(chr1_max - HG38_CHR1_LENGTH)
    tolerance = 5_000_000  # 5 Mb tolerance — AADR positions can be near the chrom end
    if dist_hg19 <= tolerance and dist_hg19 <= dist_hg38:
        return "hg19"
    if dist_hg38 <= tolerance and dist_hg38 < dist_hg19:
        return "hg38"
    raise UnsupportedAADRBuild(
        what=f"chr1 max position {chr1_max}",
        why=(
            f"matches neither hg19 ({HG19_CHR1_LENGTH}) nor hg38 ({HG38_CHR1_LENGTH}) "
            "within 5 Mb"
        ),
        fix="Pass --aadr-build hg19|hg38 to override; v0.1 expects AADR hg19 (default through v66)",
    )


__all__ = [
    "HG19_CHR1_LENGTH",
    "HG38_CHR1_LENGTH",
    "BuildOverride",
    "detect_aadr_build",
    "detect_bam_build",
    "detect_bam_format",
    "detect_bam_sample_name",
    "normalize_chrom",
    "parse_aadr_snp",
]
