"""`inspect` subcommand implementation — pure-Python AADR .snp summary.

No external dependencies (no subprocess, no network). The summary is useful for
sanity-checking an AADR panel before running `extract`, especially the
panel_guess heuristic that flags 1240k vs HO size class.
"""
from __future__ import annotations

import json
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Final

from .errors import PileupAadrError
from .format_detect import (
    classify_aadr_chrom_set,
    detect_aadr_build,
    parse_aadr_snp,
)

log = logging.getLogger(__name__)

# Strand-ambiguous (palindromic) allele pairs: a SNP whose two alleles complement each other
# is the same on forward and reverse strand. ~8% of biallelic SNPs in the 1240k panel.
_PALINDROMES: Final[frozenset[tuple[str, str]]] = frozenset(
    {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
)


def run_inspect(*, aadr_snp: Path, json_output: bool) -> int:
    """Render a structured summary of an AADR .snp file. Returns exit code (0 on success).

    Args:
        aadr_snp: path to AADR .snp file
        json_output: if True, emit JSON to stdout; else TSV

    Returns:
        0 on success; non-zero exit codes propagate via PileupAadrError raised by the parser.

    Fields:
        total_rows            int
        duplicate_rsid_count  int (always 0 — parse_aadr_snp raises on duplicates)
        build                 "hg19" | "hg38" | "unknown"
        chrom_distribution    dict[chrom_str, count]
        chrom_set             "autosomes+sex" | "autosomes_only" | "sex_only" |
                              "chrY_only" | "chrM_only" | "custom"
                              (drives whether `extract`'s autosomal coverage
                              gate applies; non-autosomal panels skip it)
        allele_distribution   dict[ref_alt_pair, count]   (e.g., "A>G": 287012)
        palindrome_count      int
        palindrome_fraction   float
        non_snp_count         int (always 0 — parse_aadr_snp raises on non-ACGT)
        morgans_present       bool
        panel_guess           "1240k" | "HO" | "unknown"
    """
    df = parse_aadr_snp(aadr_snp)  # validates duplicate-rsID + non-ACGT

    chrom_distribution = df["chrom_int"].value_counts().to_dict()
    allele_counter = Counter(zip(df["ref"], df["alt"], strict=True))
    allele_distribution = {f"{r}>{a}": int(n) for (r, a), n in allele_counter.items()}
    palindrome_count = sum(
        allele_counter[(r, a)] for (r, a) in _PALINDROMES if (r, a) in allele_counter
    )
    morgans_present = bool((df["gen_morgans"] != 0).any())

    # Panel guess: heuristic on row count
    n_rows = len(df)
    if 1_100_000 <= n_rows <= 1_300_000:
        panel_guess = "1240k"
    elif 250_000 <= n_rows <= 700_000:
        panel_guess = "HO"
    else:
        panel_guess = "unknown"

    try:
        build = detect_aadr_build(df, override="auto")
    except PileupAadrError:
        build = "unknown"

    summary: dict[str, Any] = {
        "total_rows": n_rows,
        "duplicate_rsid_count": 0,  # parse_aadr_snp would have raised
        "build": build,
        "chrom_distribution": {str(k): int(v) for k, v in chrom_distribution.items()},
        "chrom_set": classify_aadr_chrom_set(df),
        "allele_distribution": allele_distribution,
        "palindrome_count": int(palindrome_count),
        "palindrome_fraction": round(palindrome_count / n_rows, 4) if n_rows else 0.0,
        "non_snp_count": 0,  # parse_aadr_snp would have raised
        "morgans_present": morgans_present,
        "panel_guess": panel_guess,
    }

    if json_output:
        json.dump(summary, sys.stdout, indent=2, sort_keys=False)
        sys.stdout.write("\n")
    else:
        # TSV: field \t value (dict values JSON-encoded)
        sys.stdout.write("field\tvalue\n")
        for k, v in summary.items():
            if isinstance(v, dict):
                sys.stdout.write(f"{k}\t{json.dumps(v, separators=(',', ':'))}\n")
            else:
                sys.stdout.write(f"{k}\t{v}\n")
    return 0


__all__ = ["run_inspect"]
