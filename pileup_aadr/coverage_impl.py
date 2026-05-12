"""`coverage` subcommand implementation.

Wraps mosdepth to produce a per-chromosome coverage report. The output is the
diagnostic users reach for when a pre-flight `validate` PASSes but `extract`'s
coverage gate FAILs — "is the BAM low-coverage globally, or just at the AADR
1240k positions?".

mosdepth is fast (~2-3 min on 67 GB BAM) when given `--no-per-base` so it only
emits the per-chrom summary we parse here. With `--by REGIONS.bed`, it
restricts the count to the regions file (useful for AADR-only coverage).
"""
from __future__ import annotations

import json
import logging
import sys
import tempfile
from pathlib import Path
from typing import Any

from .tool_wrapper import MOSDEPTH_SPEC, ToolWrapper
from .types import CoverageCliArgs

log = logging.getLogger(__name__)


def run_coverage(args: CoverageCliArgs) -> int:
    """Per-chromosome BAM coverage report via mosdepth. Returns process exit code."""
    wrapper = ToolWrapper(MOSDEPTH_SPEC)
    with tempfile.TemporaryDirectory(prefix="pileup-aadr-coverage-") as td_str:
        td = Path(td_str)
        prefix = td / "coverage"
        mosdepth_args = [
            "--threads", str(args.threads),
            "--no-per-base",
            "--quantize", args.quantize,
        ]
        if args.regions:
            mosdepth_args.extend(["--by", str(args.regions)])
        mosdepth_args.extend([str(prefix), str(args.bam)])

        wrapper.run(
            args=mosdepth_args,
            capture_stderr_to=td / "mosdepth.stderr",
            check=True,
        )
        # mosdepth writes <prefix>.mosdepth.summary.txt with per-chrom + total stats
        summary = _parse_mosdepth_summary(
            Path(f"{prefix}.mosdepth.summary.txt")
        )

    if args.json_output:
        json.dump(summary, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        sys.stdout.write("chrom\tlength\tbases\tmean_coverage\n")
        for chrom, fields in summary["per_chrom"].items():
            sys.stdout.write(
                f"{chrom}\t{fields['length']}\t{fields['bases']}\t"
                f"{fields['mean_coverage']}\n"
            )
    return 0


def _parse_mosdepth_summary(summary_path: Path) -> dict[str, Any]:
    """Parse mosdepth's `<prefix>.mosdepth.summary.txt` into a structured dict.

    Format: TSV with columns `chrom length bases mean min max`. Mosdepth also
    emits per-region rollups (e.g., `chr1_region`) when `--by` is used; we
    keep them as-is in the per_chrom dict so the caller sees both raw +
    per-region rows. The `total` and `total_region` rollup rows are also
    preserved verbatim.
    """
    per_chrom: dict[str, dict[str, Any]] = {}
    with open(summary_path) as f:
        next(f)  # header line
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, length, bases, mean, _min, _max = parts[:6]
            per_chrom[chrom] = {
                "length": int(length),
                "bases": int(bases),
                "mean_coverage": float(mean),
            }
    return {"per_chrom": per_chrom}


__all__ = ["run_coverage"]
