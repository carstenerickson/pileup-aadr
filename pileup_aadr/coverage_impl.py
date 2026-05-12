"""`coverage` subcommand implementation.

Wraps mosdepth to produce a per-chromosome coverage report. The output is the
diagnostic users reach for when a pre-flight `validate` PASSes but `extract`'s
coverage gate FAILs — "is the BAM low-coverage globally, or just at the AADR
1240k positions?".

mosdepth is fast (~2-3 min on 67 GB BAM) when given `--no-per-base` so it only
emits the per-chrom summary we parse here. With `--by REGIONS.bed`, it
restricts the count to the regions file (useful for AADR-only coverage).

Output columns per HLD §"CLI reference > coverage":
  chrom  length  bases  mean_coverage  median_coverage
  fraction_at_>=1x  fraction_at_>=5x  fraction_at_>=10x  fraction_at_>=30x

`length`, `bases`, `mean_coverage` come from `<prefix>.mosdepth.summary.txt`;
`median_coverage` and the four `fraction_at_>=Nx` columns come from
`<prefix>.mosdepth.global.dist.txt` (cumulative depth distribution, always
written by mosdepth regardless of --quantize). The `--quantize` value is
forwarded to mosdepth so its `.quantized.bed.gz` is also written for users
who want bin-level analysis; the four threshold columns we emit are derived
from the integer-depth distribution and are stable across `--quantize`
values.
"""
from __future__ import annotations

import json
import logging
import sys
import tempfile
from pathlib import Path
from typing import Any, Final

from .tool_wrapper import MOSDEPTH_SPEC, ToolWrapper
from .types import CoverageCliArgs

log = logging.getLogger(__name__)

# Depth thresholds reported as fraction_at_>=Nx columns. Per HLD §"CLI reference
# > coverage". Stable across --quantize values (those configure mosdepth's
# .quantized.bed.gz binning, which is independent of these per-depth fractions).
_DEPTH_THRESHOLDS: Final[tuple[int, ...]] = (1, 5, 10, 30)


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
        summary = _parse_mosdepth_summary(
            Path(f"{prefix}.mosdepth.summary.txt"),
        )
        # global.dist.txt: per-chrom (depth, cumulative-fraction-at-≥-depth) rows.
        # Always written by mosdepth — independent of --quantize.
        dist = _parse_mosdepth_global_dist(
            Path(f"{prefix}.mosdepth.global.dist.txt"),
        )
        _merge_dist_into_summary(summary, dist)

    if args.json_output:
        json.dump(summary, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        cols = (
            "chrom", "length", "bases", "mean_coverage", "median_coverage",
            "fraction_at_>=1x", "fraction_at_>=5x",
            "fraction_at_>=10x", "fraction_at_>=30x",
        )
        sys.stdout.write("\t".join(cols) + "\n")
        for chrom, fields in summary["per_chrom"].items():
            sys.stdout.write(
                f"{chrom}\t{fields['length']}\t{fields['bases']}\t"
                f"{fields['mean_coverage']}\t{fields['median_coverage']}\t"
                f"{fields['fraction_at_>=1x']}\t{fields['fraction_at_>=5x']}\t"
                f"{fields['fraction_at_>=10x']}\t{fields['fraction_at_>=30x']}\n"
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


def _parse_mosdepth_global_dist(dist_path: Path) -> dict[str, list[tuple[int, float]]]:
    """Parse `<prefix>.mosdepth.global.dist.txt` into per-chrom (depth, cum-fraction).

    Format: 3-col TSV `chrom\\tdepth\\tcumulative_fraction` where rows for each
    chrom are emitted in DESCENDING depth order; cumulative_fraction is the
    fraction of bases on that chrom covered at >= depth. The `total` rollup
    is included.

    Returns:
        {chrom: [(depth, fraction), ...]} with each chrom's list sorted by
        depth descending (matches mosdepth's emission order).
    """
    per_chrom: dict[str, list[tuple[int, float]]] = {}
    if not dist_path.exists():
        # Defensive — mosdepth always writes this; missing means a broken run
        log.warning("mosdepth global.dist.txt missing at %s; quantile cols set to NaN", dist_path)
        return per_chrom
    with open(dist_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            chrom, depth_s, frac_s = parts
            per_chrom.setdefault(chrom, []).append((int(depth_s), float(frac_s)))
    return per_chrom


def _merge_dist_into_summary(
    summary: dict[str, Any],
    dist: dict[str, list[tuple[int, float]]],
) -> None:
    """Compute median_coverage + fraction_at_>=Nx columns for each chrom in summary.

    Median: largest depth where the cumulative fraction is still >= 0.5.
    fraction_at_>=Nx: cumulative fraction at depth N (lookup; default 0.0
    when N is below mosdepth's emitted minimum which means 100% — but we
    use the actual emitted value when present and default to 0.0 otherwise
    matching mosdepth's "no row → no coverage at this depth" convention).
    """
    for chrom, fields in summary["per_chrom"].items():
        rows = dist.get(chrom)
        if not rows:
            fields["median_coverage"] = float("nan")
            for n in _DEPTH_THRESHOLDS:
                fields[f"fraction_at_>={n}x"] = float("nan")
            continue
        # rows are descending by depth; the largest depth with fraction >= 0.5 is the median
        median_depth = 0
        for depth, frac in rows:
            if frac >= 0.5:
                median_depth = depth
                break
        fields["median_coverage"] = median_depth
        # Map depth → fraction for O(1) per-threshold lookup
        depth_to_frac = dict(rows)
        for n in _DEPTH_THRESHOLDS:
            fields[f"fraction_at_>={n}x"] = depth_to_frac.get(n, 0.0)


__all__ = ["run_coverage"]
