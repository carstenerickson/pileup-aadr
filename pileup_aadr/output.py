"""Run-level output writers: sidecar JSON, JSON report, per-variant TSV, stdout summary.

Most writers are thin. The actual EIGENSTRAT triplet (`<prefix>.{geno,snp,ind}`)
is written inline in `rejoin.py` for streaming reasons (1.2M-site materialization
would waste ~25 MB and the writes ARE the rejoin loop). This module owns:

- `write_pseudohaploid_sidecar` — `<prefix>.pseudohaploid.json` (consumed by pgen-samplebind)
- `write_json_report` — run-level summary JSON for the ancestry-pipeline-tool gate
- `write_per_variant_tsv` — streaming per-variant action TSV (--report-tsv)
- `write_stdout_summary` — human-readable multi-line block matching HLD §"Stdout summary"
"""
from __future__ import annotations

import dataclasses
import json
import logging
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import IO, Any

from . import __version__
from .counters import ExtractCounters

log = logging.getLogger(__name__)

JSON_REPORT_SCHEMA_VERSION = 1
PSEUDOHAPLOID_SIDECAR_SCHEMA_VERSION = 1


def write_pseudohaploid_sidecar(path: Path, sidecar: dict[str, Any]) -> None:
    """Write the `<prefix>.pseudohaploid.json` sidecar.

    Args:
        path: full sidecar path (e.g., `/data/carsten_pseudohaploid.pseudohaploid.json`).
        sidecar: dict from `rejoin.RejoinOutput.pseudohaploid_sidecar`.
    """
    if "schema_version" not in sidecar:
        sidecar = {"schema_version": PSEUDOHAPLOID_SIDECAR_SCHEMA_VERSION, **sidecar}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(sidecar, f, indent=2)
        f.write("\n")
    log.info("Wrote pseudohaploid sidecar: %s", path)


def write_json_report(
    path: Path,
    counters: ExtractCounters,
    *,
    config: dict[str, Any],
    tool_versions: dict[str, str],
    input_meta: dict[str, Any],
    output_meta: dict[str, Any],
) -> None:
    """Write the run-level summary JSON report.

    Schema-version-1 layout:
        {
          "schema_version": 1,
          "tool":   { "name", "version", "tool_versions": {...} },
          "input":  { "bam_path", "bam_format", "bam_build", ... },
          "stage_1_lift": {...},  # ExtractCounters fields (or null for fast path)
          "stage_2_transform": {...},
          "stage_3_call": {...},
          "stage_4_rejoin": {...},
          "coverage": {...},
          "gates": {...},
          "wallclock_total_seconds": float,
          "output": { "prefix", "geno_bytes", ... },
          "config": {...}
        }
    """
    payload = {
        "schema_version": JSON_REPORT_SCHEMA_VERSION,
        "tool": {
            "name": "pileup-aadr",
            "version": __version__,
            "tool_versions": tool_versions,
        },
        "input": input_meta,
        # ExtractCounters fields (stage_1_lift through wallclock_total_seconds)
        # land at the top level per the schema.
        **dataclasses.asdict(counters),
        "output": output_meta,
        "config": config,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, default=_json_default)
        f.write("\n")
    log.info("Wrote JSON report: %s", path)


def _json_default(obj: object) -> object:
    """JSON encoder for non-stdlib types (Path, set). Only fires for unusual cases."""
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, set):
        return sorted(obj)
    raise TypeError(f"Cannot serialize {type(obj).__name__}: {obj!r}")


def write_per_variant_tsv(path: Path, rows: Iterable[dict[str, Any]]) -> int:
    """Write streaming per-variant action TSV. Constant memory regardless of variant count.

    Args:
        path: output TSV path.
        rows: iterable of {"aadr_id", "chrom_hg19", "pos_hg19", "ref_hg19",
            "alt_hg19", "action"} dicts.

    Returns:
        Number of rows written.
    """
    cols = ["aadr_id", "chrom_hg19", "pos_hg19", "ref_hg19", "alt_hg19", "action"]
    n = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")
            n += 1
    log.info("Wrote per-variant TSV: %s (%d rows)", path, n)
    return n


def write_stdout_summary(
    counters: ExtractCounters,
    *,
    bam_path: Path,
    bam_format: str,
    bam_build: str,
    bam_coverage: float | None,
    aadr_path: Path,
    aadr_total: int,
    ref_fasta: Path,
    chain_path: Path,
    output_prefix: Path,
    output_bytes: dict[str, int],
    sex: str = "U",
) -> None:
    """Render the human-readable multi-line summary block to stdout.

    Mirrors HLD §"Stdout summary" example layout. Suppressed by --quiet
    (orchestrator skips this call when quiet=True).
    """
    out = sys.stdout

    cov_str = f", {bam_coverage:.1f}x mean coverage" if bam_coverage is not None else ""
    out.write(f"pileup-aadr extract {bam_path} {aadr_path} -o {output_prefix}\n\n")
    out.write("Input:\n")
    out.write(f"  BAM:        {bam_path} ({bam_format}, {bam_build}{cov_str})\n")
    out.write(f"  AADR .snp:  {aadr_path} ({aadr_total:,} variants)\n")
    out.write(f"  Reference:  {ref_fasta}\n")
    out.write(f"  Chain:      {chain_path}\n\n")

    if counters.stage_1_lift is not None:
        s1 = counters.stage_1_lift
        out.write(
            f"Stage 1 - lift AADR sites hg19 -> hg38 (Picard LiftoverVcf "
            f"RECOVER_SWAPPED_REF_ALT) [{s1.wallclock_seconds:.1f}s]\n"
        )
        out.write(f"  Input variants:                 {aadr_total:>10,}\n")
        out.write(
            f"  Pre-lift filter (palindromes):  "
            f"{-s1.input_filters['palindrome_drops']:>10,}\n"
        )
        out.write(
            f"  Pre-lift filter (non-SNP):      "
            f"{-s1.input_filters['non_snp_drops']:>10,}\n"
        )
        for reason, count in s1.rejected_by_reason.items():
            if count > 0:
                out.write(f"  Picard rejected ({reason}):     {-count:>10,}\n")
        out.write(
            f"  Picard 'SwappedAlleles' INFO marker: {s1.swapped_alleles_count:,} "
            f"(recovered, not rejected)\n"
        )
        gate = "WARN" if s1.liftover_yield_warning else "PASS"
        out.write(
            f"  Lifted: {s1.lifted_sites:,} sites "
            f"(yield: {s1.liftover_yield_pct:.2f}% - gate {gate})\n\n"
        )

    if counters.stage_2_transform is not None:
        s2 = counters.stage_2_transform
        out.write(
            f"Stage 2 - transform to pileupCaller .snp [{s2.wallclock_seconds:.1f}s]\n"
        )
        out.write(f"  Alt-contig filter (chrN_*_alt): {-s2.alt_contig_drops:>10,}\n")
        out.write(f"  Output: {s2.output_sites:,} sites\n\n")

    s3 = counters.stage_3_call
    pc = s3.pileupcaller_summary
    minutes, seconds = divmod(s3.wallclock_seconds, 60)
    s3_time = f"{int(minutes)}m {int(seconds)}s"
    out.write(
        f"Stage 3 - samtools mpileup + pileupCaller --randomDiploid [{s3_time}]\n"
    )
    out.write("  pileupCaller stderr summary stats (parsed for JSON report):\n")
    out.write(f"    TotalSites:        {pc.total_sites:,}\n")
    out.write(f"    NonMissingCalls:   {pc.non_missing_calls:,}\n")
    out.write(f"    avgRawReads:       {pc.avg_raw_reads:.1f}\n")
    out.write(f"    avgSampledFrom:    {pc.avg_sampled_from:.1f}\n\n")

    if counters.stage_4_rejoin is not None:
        s4 = counters.stage_4_rejoin
        out.write(
            f"Stage 4 - rejoin hg19 coordinates by rsID + invert dosage at "
            f"SwappedAlleles [{s4.wallclock_seconds:.1f}s]\n"
        )
        out.write(f"  rsID matched in AADR .snp:        {s4.rsid_matched:,}\n")
        if s4.rsid_matched > 0:
            swap_pct = 100.0 * s4.ref_alt_swap_count / s4.rsid_matched
            out.write(
                f"  REF/ALT swapped (dosage inverted): {s4.ref_alt_swap_count:,} "
                f"({swap_pct:.2f}%)\n"
            )
        out.write(
            f"  Allele mismatch (dropped):         {s4.allele_mismatch_drops:,}  "
            f"(defensive sanity check; ~0 expected)\n"
        )
        out.write(f"  Output: {s4.output_variants:,} variants x 1 sample\n\n")

    cov = counters.coverage
    out.write("Coverage report (post-rejoin):\n")
    out.write(
        f"  Total non-missing autosomal calls: {cov.non_missing_autosomal_calls:,} "
        f"of {aadr_total:,}  ({cov.coverage_fraction * 100:.1f}%)\n"
    )
    out.write("  Per-chromosome (all 22 autosomes; gated count):\n")
    _write_chrom_grid(out, cov.per_chrom_call_count)
    out.write("  Sex chromosomes + mtDNA (informational; not in gated count):\n")
    chrx = cov.per_chrom_call_count.get("chrX", 0)
    chry = cov.per_chrom_call_count.get("chrY", 0)
    chrm = cov.per_chrom_call_count.get("chrM", 0)
    sex_note = "(PSEUDOHAPLOID male/U sex)" if sex in ("M", "U") and chrx == 0 else ""
    out.write(f"    chrX:   {chrx:>8,}   {sex_note}\n")
    out.write(f"    chrY:   {chry:>8,}\n")
    out.write(f"    chrM:   {chrm:>8,}   (not in AADR 1240k panel)\n")

    cov_gate = counters.gates.get("coverage", "PASS")
    out.write(f"  Coverage gate: >=800k threshold -> {cov_gate}\n\n")

    total_bytes = sum(output_bytes.values())
    out.write(
        f"Wrote {output_prefix}.{{geno,snp,ind,pseudohaploid.json}} "
        f"({_human_bytes(total_bytes)} total)\n"
    )
    minutes, seconds = divmod(counters.wallclock_total_seconds, 60)
    out.write(f"Done in {int(minutes)}m {int(seconds)}s.\n")
    out.flush()


def _write_chrom_grid(out: IO[str], per_chrom: dict[str, int]) -> None:
    """Write 22 autosomes in a 6-row x 4-col grid for readability."""
    rows = 6
    autosomes = [f"chr{i}" for i in range(1, 23)]
    cols = [autosomes[i:i + rows] for i in range(0, len(autosomes), rows)]
    max_chroms_per_row = max(len(c) for c in cols)
    for r in range(max_chroms_per_row):
        line_parts = []
        for c in cols:
            if r < len(c):
                chrom = c[r]
                count = per_chrom.get(chrom, 0)
                line_parts.append(f"{chrom}: {count:>7,}")
        out.write("    " + "   ".join(line_parts) + "\n")


def _human_bytes(n: int) -> str:
    """Format byte count as human-readable string (B/KB/MB/GB/TB)."""
    size = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024:
            return f"{int(size)} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


__all__ = [
    "JSON_REPORT_SCHEMA_VERSION",
    "PSEUDOHAPLOID_SIDECAR_SCHEMA_VERSION",
    "write_json_report",
    "write_per_variant_tsv",
    "write_pseudohaploid_sidecar",
    "write_stdout_summary",
]
