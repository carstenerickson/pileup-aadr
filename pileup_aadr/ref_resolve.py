"""Reference FASTA resolution + build-match verification.

Per HLD §"Chain & reference dependencies" + LLD §15 (extract_orch helpers):
  - Resolution order: --ref-fasta → $PILEUP_AADR_REF_DIR/<build>.fa → BAM @PG extraction
  - Build verification: read .fai chr1 length; compare against BAM-detected build
    (±1 Mb tolerance) — catches wrong-build FASTA before Picard wastes ~10s loading
    + ~30s rejecting most sites

This module exists as its own file (not folded into extract_orch.py) so validate.py
can use it for pre-flight checks without depending on the full orchestrator.
"""
from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Final, Literal

import pysam

from .errors import ReferenceFastaBuildMismatchError, ReferenceFastaNotFound
from .format_detect import HG19_CHR1_LENGTH, HG38_CHR1_LENGTH

log = logging.getLogger(__name__)

_BUILD_TOLERANCE_BP: Final[int] = 1_000_000  # ±1 Mb tolerance for patch-level differences

_CHR1_LENGTH_BY_BUILD: Final[dict[str, int]] = {
    "hg19": HG19_CHR1_LENGTH,
    "hg38": HG38_CHR1_LENGTH,
}


def resolve_ref_fasta(
    cli_ref: Path | None,
    bam: Path,
    bam_build: Literal["hg19", "hg38"],
) -> Path:
    """3-tier resolution of the target FASTA per HLD §"Chain & reference dependencies".

    Resolution order:
      1. --ref-fasta PATH (explicit override)
      2. $PILEUP_AADR_REF_DIR/<build>.fa env var
      3. Auto-detect from BAM @PG headers (samtools/bwa/dragen all record reference path)

    Every resolved path is passed through `verify_fasta_matches_bam_build` before
    return — explicit --ref-fasta is NOT exempt (catches wrong-build FASTA before
    Stage 1 wastes time).

    Args:
        cli_ref: --ref-fasta PATH if user passed it
        bam: aligned BAM/CRAM (only used for @PG fallback extraction)
        bam_build: "hg19" or "hg38" — the BAM's detected build, drives the
            $PILEUP_AADR_REF_DIR/<build>.fa lookup + the FASTA verification

    Returns:
        Path to a verified target FASTA matching `bam_build`.

    Raises:
        ReferenceFastaNotFound: no resolution path produced a readable FASTA.
        ReferenceFastaBuildMismatchError: the resolved FASTA's chr1 length disagrees
            with `bam_build`'s canonical value beyond ±1 Mb tolerance.
    """
    if cli_ref is not None:
        verify_fasta_matches_bam_build(cli_ref, bam_build)
        return cli_ref

    env_dir = os.environ.get("PILEUP_AADR_REF_DIR")
    if env_dir:
        candidate = Path(env_dir) / f"{bam_build}.fa"
        if candidate.exists():
            verify_fasta_matches_bam_build(candidate, bam_build)
            return candidate

    extracted = _extract_ref_from_bam_pg(bam)
    verify_fasta_matches_bam_build(extracted, bam_build)
    return extracted


def _extract_ref_from_bam_pg(bam: Path) -> Path:
    """Parse BAM @PG header lines for the alignment's reference FASTA path.

    samtools, bwa, dragen, and others record the reference path in @PG fields
    (typically `CL:` for the command line or `DS:` for description). Search for
    any path-looking value ending in `.fa` or `.fasta` and return the first one
    that exists on disk.

    Raises:
        ReferenceFastaNotFound: no readable .fa path found in any @PG record.
    """
    fasta_re = re.compile(r"(/[^\s]+\.(?:fa|fasta))")
    with pysam.AlignmentFile(str(bam), "r", check_sq=False) as bf:
        pgs = bf.header.to_dict().get("PG", [])

    candidates: list[str] = []
    for pg in pgs:
        for value in pg.values():
            if isinstance(value, str):
                candidates.extend(m.group(1) for m in fasta_re.finditer(value))

    for cand in candidates:
        if Path(cand).exists():
            log.debug("Resolved ref FASTA via BAM @PG: %s", cand)
            return Path(cand)

    raise ReferenceFastaNotFound(
        what=str(bam),
        why=(
            f"no readable .fa path found in BAM @PG records "
            f"({len(candidates)} candidates inspected)"
        ),
        fix=(
            "Pass --ref-fasta PATH explicitly, or set $PILEUP_AADR_REF_DIR. "
            f"@PG candidates seen (none readable): {candidates[:3]!r}"
        ),
    )


def verify_fasta_matches_bam_build(fasta: Path, bam_build: Literal["hg19", "hg38"]) -> None:
    """Read the FASTA's .fai chr1 length; raise if it disagrees with `bam_build`.

    Cheap pre-flight (~1 ms — just reads the .fai) that catches the wrong-build FASTA
    case BEFORE Picard wastes ~10 sec loading the reference + ~30 sec rejecting
    most sites with MismatchedRefAllele.

    If the .fai is missing (uncommon — most BAMs come alongside a faidx'd FASTA),
    we skip the check rather than raising; Picard will discover the issue on its own
    and emit its own diagnostic.

    Raises:
        ReferenceFastaBuildMismatchError: chr1 length disagrees beyond ±1 Mb tolerance.
    """
    fai_path = Path(f"{fasta}.fai")
    if not fai_path.exists():
        log.debug("Skipping FASTA build verification: %s missing", fai_path)
        return

    chr1_length: int | None = None
    with open(fai_path) as f:
        for line in f:
            parts = line.split("\t")
            if not parts:
                continue
            name = parts[0]
            if name in ("chr1", "1") and len(parts) >= 2:
                chr1_length = int(parts[1])
                break

    if chr1_length is None:
        log.debug("Skipping FASTA build verification: no chr1/1 record in %s", fai_path)
        return

    # Closest-match: hg19 chr1 (249,250,621) and hg38 chr1 (248,956,422) are only
    # 294 KB apart, so the ±1 Mb tolerance window around either build INCLUDES the
    # other build's exact chr1 length. We pick the build whose canonical chr1 is
    # closest to the FASTA's, then verify that closest match also fits the
    # patch-level tolerance window.
    distances = {b: abs(chr1_length - exp) for b, exp in _CHR1_LENGTH_BY_BUILD.items()}
    closest_build, closest_dist = min(distances.items(), key=lambda kv: kv[1])

    if closest_build != bam_build or closest_dist > _BUILD_TOLERANCE_BP:
        # Either the FASTA matches a different build, or it matches no known build
        # (e.g., T2T-CHM13). Diagnose with whichever's closer.
        raise ReferenceFastaBuildMismatchError(
            what=str(fasta),
            why=(
                f"FASTA chr1 length {chr1_length} matches {closest_build} "
                f"(distance {closest_dist}) but BAM is {bam_build}"
            ),
            fix=(
                f"Pass --ref-fasta pointing at a {bam_build} FASTA; or override with "
                f"--bam-build {closest_build} if the BAM is actually {closest_build}-aligned"
            ),
        )
    log.debug(
        "FASTA build verified: %s matches BAM build %s (chr1 length %d)",
        fasta,
        bam_build,
        chr1_length,
    )


__all__ = [
    "resolve_ref_fasta",
    "verify_fasta_matches_bam_build",
]
