"""Target FASTA `.dict` resolution + auto-generation.

Per HLD §"Chain & reference dependencies > Target FASTA sequence dictionary":

  1. If `<fasta>.dict` exists → use it
  2. Else if `<fasta-stem>.dict` exists → use it
  3. Else generate via `picard CreateSequenceDictionary R=<fasta> O=<dict>` and
     cache. Cache target: alongside the FASTA if its parent dir is writable;
     otherwise `~/.cache/pileup-aadr/dicts/<sha-of-abs-path>.dict`.

One-time cost on hg38: ~23 sec (verified empirically v2.1); subsequent runs
reuse the cached dict.

Picard requires a `.dict` for any LiftoverVcf call against a fresh FASTA;
without one, it errors at startup with an opaque message about "missing
sequence dictionary". This module catches that case at extract pre-flight
+ validate time, so users get the auto-generated path and not the failure.
"""
from __future__ import annotations

import hashlib
import logging
import os
from pathlib import Path

from .tool_wrapper import PICARD_SPEC, ToolWrapper

log = logging.getLogger(__name__)

_USER_CACHE_DIR_DEFAULT = Path.home() / ".cache" / "pileup-aadr" / "dicts"


def find_existing_dict(fasta: Path) -> Path | None:
    """Return the `.dict` path alongside `fasta`, or None if neither candidate exists.

    Picard accepts both `<fasta>.dict` (e.g., `hg38.fa.dict`) and `<fasta-stem>.dict`
    (e.g., `hg38.dict`). The latter is the convention shipped with most reference
    bundles (GATK resource bundle, bcbio-nextgen, etc.); the former is what
    Picard's CreateSequenceDictionary writes by default.
    """
    long_form = Path(f"{fasta}.dict")  # <fasta>.dict — e.g., hg38.fa.dict
    if long_form.exists():
        return long_form
    short_form = fasta.with_suffix(".dict")  # <fasta-stem>.dict — e.g., hg38.dict
    if short_form.exists():
        return short_form
    return None


def find_or_user_cache_dict_path(fasta: Path, *, cache_dir: Path | None = None) -> Path:
    """Decide where to write a NEW `.dict` for `fasta`.

    Prefer alongside the FASTA when its parent dir is writable (most users have
    write access to their own data dir). Fall back to a user-cache location
    keyed by the SHA-256 of the FASTA's absolute path so two FASTAs with the
    same basename don't collide.

    Args:
        fasta: target FASTA path
        cache_dir: override for the user-cache location (default
            `~/.cache/pileup-aadr/dicts/`); used by tests for isolation.
    """
    fasta_abs = fasta.resolve()
    long_form = Path(f"{fasta_abs}.dict")
    if os.access(fasta_abs.parent, os.W_OK):
        return long_form
    cache = cache_dir if cache_dir is not None else _USER_CACHE_DIR_DEFAULT
    cache.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(str(fasta_abs).encode()).hexdigest()[:16]
    return cache / f"{fasta_abs.name}.{digest}.dict"


def ensure_target_fasta_dict(
    fasta: Path,
    *,
    cache_dir: Path | None = None,
) -> Path:
    """Resolve or generate the `.dict` sidecar required by Picard.

    Returns:
        Path to a `.dict` file that exists on disk and matches `fasta`.

    Raises:
        ToolSubprocessError: Picard's CreateSequenceDictionary subprocess failed.
        PicardNotFoundError: Picard JAR unresolvable (and we needed to generate).
    """
    existing = find_existing_dict(fasta)
    if existing is not None:
        log.debug("Found existing .dict for %s at %s", fasta, existing)
        return existing

    out = find_or_user_cache_dict_path(fasta, cache_dir=cache_dir)
    log.info(
        "Generating sequence dictionary for %s -> %s (one-time ~23s on hg38)",
        fasta, out,
    )
    out.parent.mkdir(parents=True, exist_ok=True)

    wrapper = ToolWrapper(PICARD_SPEC)
    stderr_path = out.parent / f"{out.name}.picard.stderr"
    wrapper.run(
        args=[
            "CreateSequenceDictionary",
            "--REFERENCE", str(fasta),
            "--OUTPUT", str(out),
        ],
        capture_stderr_to=stderr_path,
        check=True,
    )
    log.info("Wrote sequence dictionary: %s", out)
    return out


__all__ = [
    "ensure_target_fasta_dict",
    "find_existing_dict",
    "find_or_user_cache_dict_path",
]
