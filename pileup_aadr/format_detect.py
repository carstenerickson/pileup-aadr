"""Pre-stage detection of input formats + builds.

Cheapest, fastest module — pure metadata reads, no genotype I/O. Called once at
orchestrator startup before any expensive work begins. Failures surface immediately
with clear errors per HLD §"Error class taxonomy".
"""
from __future__ import annotations

import hashlib
import logging
import os
import re
from pathlib import Path
from typing import Final, Literal

import pandas as pd
import pysam

from .chrom_lengths import CHROM_ORDER
from .errors import (
    AADRDuplicateRsidError,
    AADRParseError,
    BAMParseError,
    BAMSampleNameAmbiguous,
    UnsupportedAADRBuild,
    UnsupportedReferenceBuild,
)

log = logging.getLogger(__name__)

# Bump whenever the parse output schema changes (column set, dtypes, sort order).
# Stale cache entries with a different version are automatically ignored (different filename).
PARSE_SCHEMA_VERSION: Final[int] = 1

# chr-prefixed chrom → sort position; used for pre-sorting parse output.
_CHROM_SORT_MAP: Final[dict[str, int]] = {c: i for i, c in enumerate(CHROM_ORDER)}

BuildOverride = Literal["hg19", "hg38", "auto"]

# Canonical chromosome lengths per assembly (hardcoded; chrom lengths are byte-stable).
# chr1 is the primary anchor; chr20 is the fallback when chr1 is absent (e.g.,
# chrY-only BAMs in haplogroup workflows). chr20 has a wider hg19/hg38 gap
# (~1.4 Mb) than chr1's 294 KB, so the per-anchor tolerance lookup picks the
# correct build cleanly via closest-match.
HG19_CHR1_LENGTH: Final[int] = 249_250_621
HG38_CHR1_LENGTH: Final[int] = 248_956_422
HG19_CHR20_LENGTH: Final[int] = 63_025_520
HG38_CHR20_LENGTH: Final[int] = 64_444_167
_BUILD_TOLERANCE_BP: Final[int] = 1_000_000  # ±1Mb tolerance for patch-level differences

# Per-anchor build lookup: tried in order until one anchor's chrom name is
# present. Each entry: (anchor_name, {build: canonical_length}).
_BUILD_ANCHORS: Final[tuple[tuple[tuple[str, ...], dict[str, int]], ...]] = (
    (("chr1", "1"), {"hg19": HG19_CHR1_LENGTH, "hg38": HG38_CHR1_LENGTH}),
    (("chr20", "20"), {"hg19": HG19_CHR20_LENGTH, "hg38": HG38_CHR20_LENGTH}),
)


def _closest_build_for(
    observed_length: int, anchors: dict[str, int],
) -> tuple[str, int]:
    """Return (build, distance) for the build whose canonical length is
    closest to `observed_length`. Caller checks distance vs tolerance."""
    distances = {b: abs(observed_length - exp) for b, exp in anchors.items()}
    return min(distances.items(), key=lambda kv: kv[1])

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
        sn_to_ln = {sq["SN"]: int(sq["LN"]) for sq in sq_records}

    # Try anchors in order — chr1 first, chr20 fallback for chrY-only BAMs etc.
    # Closest-match per anchor: hg19 chr1 (249,250,621) and hg38 chr1
    # (248,956,422) are only 294 KB apart — well INSIDE the ±1 Mb tolerance,
    # so first-match-wins logic returned "hg19" for every hg38 BAM (#1).
    # Closest-match disambiguates cleanly.
    last_anchor_seen: str | None = None
    for anchor_names, anchors in _BUILD_ANCHORS:
        observed = next(
            (sn_to_ln[name] for name in anchor_names if name in sn_to_ln), None,
        )
        if observed is None:
            continue
        last_anchor_seen = anchor_names[0]
        closest_build, closest_dist = _closest_build_for(observed, anchors)
        if closest_dist <= _BUILD_TOLERANCE_BP:
            return closest_build  # type: ignore[return-value]
        # The anchor is present but its length doesn't match either build —
        # likely a non-hg19/hg38 assembly (T2T-CHM13). No point trying further
        # anchors; they'd produce the same diagnosis.
        raise UnsupportedReferenceBuild(
            what=str(bam_path),
            why=(
                f"BAM @SQ {anchor_names[0]} length {observed} matches neither "
                f"hg19 ({anchors['hg19']}) nor hg38 ({anchors['hg38']}) within ±1 Mb"
            ),
            fix=(
                f"Pass --bam-build hg19|hg38 to override (if the BAM is "
                f"hg19/hg38-compatible despite differing {anchor_names[0]} length, "
                f"e.g., T2T-CHM13 fragments)"
            ),
        )

    # Fell through both anchors — no chr1 OR chr20 in @SQ. Diagnostic names
    # both anchors so users with esoteric chrom-only BAMs know what we tried.
    if last_anchor_seen is None:
        raise UnsupportedReferenceBuild(
            what=str(bam_path),
            why=(
                "BAM @SQ has neither chr1/1 nor chr20/20 record; cannot "
                "determine build"
            ),
            fix="Pass --bam-build hg19|hg38 to override detection",
        )
    # Unreachable — `last_anchor_seen` would have triggered the per-anchor
    # raise above. Placeholder for the type-checker.
    raise UnsupportedReferenceBuild(  # pragma: no cover
        what=str(bam_path),
        why="internal error: anchor detection fell through",
        fix="Pass --bam-build hg19|hg38 to override",
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


def _aadr_cache_path(sha256: str) -> Path:
    """XDG-respecting path for the schema-versioned feather cache of a parsed .snp."""
    cache_home = Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache")))
    return cache_home / "pileup-aadr" / "snp" / f"{sha256}.v{PARSE_SCHEMA_VERSION}.feather"


def parse_aadr_snp(aadr_snp_path: Path) -> pd.DataFrame:
    """Load AADR .snp into a DataFrame indexed by rsid, sorted by (CHROM_ORDER, pos_bp).

    Results are content-addressed cached as feather at
    $XDG_CACHE_HOME/pileup-aadr/snp/<sha256>.v{PARSE_SCHEMA_VERSION}.feather.
    Cache hits skip the parse + sort entirely (~3s → ~100ms on the 1.1M-site panel).

    Format (verified empirically v2.1):
        6 columns whitespace-separated (multi-space-padded for column alignment, NOT
        tab-separated). No header. Lines starting with '#' tolerated as comments.

    Args:
        aadr_snp_path: AADR .snp file

    Returns:
        DataFrame with 6 typed columns + rsid index, sorted by (CHROM_ORDER, pos_bp):
            rsid: str (index)
            chrom_int: str (kept as string to handle "1"-"22","23","24","90","91")
            gen_morgans: float
            pos_bp: int
            ref: str (single ACGT char)
            alt: str (single ACGT char)
        Rows with chrom_int not in CHROM_ORDER sort to the end (stable).

    Raises:
        AADRParseError: row has wrong column count, non-ACGT alleles, or unparseable position
        AADRDuplicateRsidError: duplicate rsid in file (Stage 4 join requires unique IDs)
    """
    sha256 = hashlib.sha256(aadr_snp_path.read_bytes()).hexdigest()
    cache_path = _aadr_cache_path(sha256)
    if cache_path.exists():
        log.info("AADR .snp cache hit (v%d): %s", PARSE_SCHEMA_VERSION, cache_path)
        return pd.read_feather(cache_path).set_index("rsid")

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

    # Sort by CHROM_ORDER x pos_bp so downstream consumers (sites_vcf, _write_aadr_native)
    # get ordered output without re-sorting. chrom_int values not in CHROM_ORDER sort last.
    chrom_rank = df["chrom_int"].map(normalize_chrom).map(_CHROM_SORT_MAP)
    chrom_rank_filled = chrom_rank.fillna(len(CHROM_ORDER))
    df = df.assign(_ck=chrom_rank_filled).sort_values(["_ck", "pos_bp"]).drop(columns="_ck")

    log.info(
        "Loaded AADR .snp: %d rows, %d unique chromosomes",
        len(df), df["chrom_int"].nunique(),
    )

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.reset_index().to_feather(cache_path)
        log.info("AADR .snp cached (v%d): %s", PARSE_SCHEMA_VERSION, cache_path)
    except Exception as exc:
        log.debug("AADR parse cache write skipped: %s", exc)

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

    # Per-anchor lookup: try chr1 (numeric "1" in AADR encoding), fall back to
    # chr20 ("20") for the chrY-only / chr-only AADR slice case. AADR positions
    # are typically near the chrom end, so the per-anchor tolerance is wider
    # (5 Mb) than the BAM @SQ length check.
    aadr_anchors = (
        ("1", {"hg19": HG19_CHR1_LENGTH, "hg38": HG38_CHR1_LENGTH}),
        ("20", {"hg19": HG19_CHR20_LENGTH, "hg38": HG38_CHR20_LENGTH}),
    )
    aadr_tolerance = 5_000_000

    last_anchor_seen: str | None = None
    for anchor_chrom, anchors in aadr_anchors:
        rows = aadr_df[aadr_df["chrom_int"] == anchor_chrom]
        if rows.empty:
            continue
        last_anchor_seen = anchor_chrom
        observed_max = int(rows["pos_bp"].max())
        closest_build, closest_dist = _closest_build_for(observed_max, anchors)
        if closest_dist <= aadr_tolerance:
            return closest_build  # type: ignore[return-value]
        # Present-but-unrecognized — same logic as detect_bam_build's per-anchor raise
        raise UnsupportedAADRBuild(
            what=f"chr{anchor_chrom} max position {observed_max}",
            why=(
                f"matches neither hg19 ({anchors['hg19']}) nor hg38 "
                f"({anchors['hg38']}) within {aadr_tolerance // 1_000_000} Mb"
            ),
            fix="Pass --aadr-build hg19|hg38 to override detection",
        )

    if last_anchor_seen is None:
        raise UnsupportedAADRBuild(
            what="(AADR DataFrame)",
            why=(
                "no chr1 (chrom_int='1') OR chr20 (chrom_int='20') rows found; "
                "cannot determine build"
            ),
            fix="Pass --aadr-build hg19|hg38 to override detection",
        )

    # Unreachable — `last_anchor_seen` would have triggered the per-anchor raise
    raise UnsupportedAADRBuild(  # pragma: no cover
        what="(AADR DataFrame)",
        why="internal error: anchor detection fell through",
        fix="Pass --aadr-build hg19|hg38 to override",
    )


AADR_CHROM_AUTOSOMES: Final[frozenset[str]] = frozenset(str(i) for i in range(1, 23))
AADR_CHROM_SEX: Final[frozenset[str]] = frozenset({"23", "24"})  # X, Y in AADR encoding
AADR_CHROM_MT: Final[frozenset[str]] = frozenset({"90"})  # mtDNA in AADR encoding


def classify_aadr_chrom_set(aadr_df: pd.DataFrame) -> str:
    """Classify the chrom set of an AADR `.snp` panel.

    Returns one of:
        "autosomes+sex"   typical 1240k / HO panels (autosomes + X + Y, optional MT)
        "autosomes_only"  full 1240k autosomal subset
        "sex_only"        chrX + chrY only (e.g., sex-chrom-only AADR slice)
        "chrY_only"       Y-chromosome haplogroup workflows
        "chrM_only"       mtDNA-only ancient workflows
        "custom"          anything else (e.g., chr22-only test slice; unusual panels)

    The classification drives whether `extract`'s autosomal coverage gate
    applies. Non-autosomal panels (chrY-only, chrM-only) skip the gate
    cleanly rather than failing with a misleading "below --min-coverage"
    message — the autosomal threshold is meaningless for haplogroup or
    mtDNA workflows. v0.2 enhancement (HLD §"Out-of-scope reachable
    extensions").
    """
    chroms_present = set(aadr_df["chrom_int"].astype(str).unique())
    has_autosome = bool(chroms_present & AADR_CHROM_AUTOSOMES)
    has_sex = bool(chroms_present & AADR_CHROM_SEX)
    has_mt = bool(chroms_present & AADR_CHROM_MT)

    if has_autosome and has_sex:
        return "autosomes+sex"
    if has_autosome and not has_sex and not has_mt:
        # All autosomes present is "autosomes_only"; partial autosomes (e.g.,
        # a chr22-only slice) is "custom" — distinguishes 1240k-autosomal-subset
        # from a one-chrom test panel.
        if chroms_present == AADR_CHROM_AUTOSOMES:
            return "autosomes_only"
        return "custom"
    if chroms_present == AADR_CHROM_SEX:
        return "sex_only"
    if chroms_present == {"24"}:  # only chrY
        return "chrY_only"
    if chroms_present == {"90"}:  # only chrM
        return "chrM_only"
    return "custom"


__all__ = [
    "AADR_CHROM_AUTOSOMES",
    "AADR_CHROM_MT",
    "AADR_CHROM_SEX",
    "HG19_CHR1_LENGTH",
    "HG19_CHR20_LENGTH",
    "HG38_CHR1_LENGTH",
    "HG38_CHR20_LENGTH",
    "PARSE_SCHEMA_VERSION",
    "BuildOverride",
    "classify_aadr_chrom_set",
    "detect_aadr_build",
    "detect_bam_build",
    "detect_bam_format",
    "detect_bam_sample_name",
    "normalize_chrom",
    "parse_aadr_snp",
]
