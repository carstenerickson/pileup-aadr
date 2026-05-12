"""`extract` subcommand orchestrator: sequences Stages 1-4 + gates + writers.

Sequence (HLD §"Module orchestration"):

  1. Pre-flight: format detection, dependency resolution, output-prefix conflict check
  2. Tool version-detection (samtools, pileupCaller, picard+java if lift path)
  3. Branch on (bam_build, aadr_build) for fast-path vs full-path
  4. Stages: 1 -> yield-gate -> 2 -> 3 -> 4 (or just 3 for fast path)
  5. Coverage gate evaluation
  6. Output writers: triplet (already inline in Stage 4), sidecar, JSON, TSV, stdout
  7. Tempdir + lock cleanup (handled by context managers)

The CLI subcommand in `extract_cmd.py` is a thin click wrapper that constructs
`ExtractCliArgs` and delegates to `run_extract`.
"""
from __future__ import annotations

import dataclasses
import hashlib
import logging
import os
import time
from pathlib import Path
from typing import Any

from . import format_detect, lift, output, pileup_call, rejoin, sites_vcf, transform
from .concurrency import output_lock, tempdir, warn_if_networked_fs
from .counters import ExtractCounters
from .errors import (
    CoverageGateFailure,
    OutputExistsError,
    PileupAadrError,
)
from .lift import chain_file_path
from .ref_resolve import resolve_ref_fasta
from .tool_wrapper import (
    JAVA_SPEC,
    PICARD_SPEC,
    PILEUPCALLER_SPEC,
    SAMTOOLS_SPEC,
    ToolWrapper,
)
from .types import ExtractCliArgs

log = logging.getLogger(__name__)


def run_extract(args: ExtractCliArgs) -> int:
    """Orchestrator for the `extract` subcommand. Returns process exit code.

    Returns:
        0 on success/warning; otherwise the `exit_code` from the raised
        `PileupAadrError`.

    Raises:
        PileupAadrError: any pileup-aadr error propagates out for `cli.py`'s
            outer formatter. All other exceptions also propagate (formatted
            as crashes by `cli.py`).
    """
    t_start = time.perf_counter()

    warn_if_networked_fs(args.output_prefix)

    bam_format = format_detect.detect_bam_format(args.bam)
    bam_build = format_detect.detect_bam_build(args.bam, override=args.bam_build)
    sample_name = format_detect.detect_bam_sample_name(
        args.bam, explicit=args.sample_name
    )
    pop_name = args.pop_name or sample_name

    aadr_df = format_detect.parse_aadr_snp(args.aadr_snp)
    aadr_build = format_detect.detect_aadr_build(aadr_df, override=args.aadr_build)

    chain = chain_file_path(
        cli_chain=args.chain_path,
        env_chain_dir=(
            Path(os.environ["PILEUP_AADR_CHAIN_DIR"])
            if "PILEUP_AADR_CHAIN_DIR" in os.environ
            else None
        ),
        strict_sha=args.strict_chain_sha,
        insecure=args.insecure_chain,
    )
    ref_fasta = resolve_ref_fasta(
        cli_ref=args.ref_fasta, bam=args.bam, bam_build=bam_build,
    )

    no_lift = bam_build == aadr_build
    log.info(
        "Pre-flight OK: BAM=%s (%s, %s), AADR=%s (%s), %s",
        args.bam, bam_format, bam_build, args.aadr_snp, aadr_build,
        "no-lift fast path" if no_lift else "full lift path",
    )

    tool_versions: dict[str, str] = {}
    tool_versions["samtools"] = ToolWrapper(SAMTOOLS_SPEC).version()
    tool_versions["pileupCaller"] = ToolWrapper(PILEUPCALLER_SPEC).version()
    if not no_lift:
        tool_versions["java"] = ToolWrapper(JAVA_SPEC).version()
        tool_versions["picard"] = ToolWrapper(PICARD_SPEC).version()

    try:
        with (
            output_lock(args.output_prefix),
            tempdir(
                base=args.tempdir,
                keep_always=args.keep_tempdir,
                clean_on_crash=args.clean_tempdir_on_crash,
            ) as td,
        ):
            # H6 fix: output-prefix existence check happens INSIDE the locked
            # region so no other process can write between check and our writes.
            _check_output_prefix(args.output_prefix, args.overwrite)
            (lift_dir := td / "lift").mkdir(exist_ok=True)
            (transform_dir := td / "transform").mkdir(exist_ok=True)
            (call_dir := td / "call").mkdir(exist_ok=True)

            counters, rejoin_out = _run_stages(
                args=args, aadr_df=aadr_df, aadr_build=aadr_build,
                chain=chain, ref_fasta=ref_fasta,
                sample_name=sample_name, pop_name=pop_name, no_lift=no_lift,
                td_lift=lift_dir, td_transform=transform_dir, td_call=call_dir,
            )

        _evaluate_coverage_gate(counters, args)
        counters = _finalize_counters(counters, time.perf_counter() - t_start)

        _write_outputs(
            args=args, counters=counters, rejoin_out=rejoin_out,
            bam_format=bam_format, bam_build=bam_build,
            aadr_total=len(aadr_df), aadr_build=aadr_build,
            ref_fasta=ref_fasta, chain=chain,
            tool_versions=tool_versions,
        )

        log.info("Extract complete in %.1fs", counters.wallclock_total_seconds)
        return 0

    except PileupAadrError as e:
        log.error("Extract failed: %s: %s", type(e).__name__, e.why)
        raise


def _check_output_prefix(prefix: Path, overwrite: bool) -> None:
    """Refuse to clobber existing outputs unless --overwrite."""
    candidates = [
        Path(f"{prefix}.geno"), Path(f"{prefix}.snp"),
        Path(f"{prefix}.ind"), Path(f"{prefix}.pseudohaploid.json"),
    ]
    existing = [p for p in candidates if p.exists()]
    if existing and not overwrite:
        raise OutputExistsError(
            what=", ".join(str(p) for p in existing),
            why=f"{len(existing)} output file(s) already exist at prefix {prefix}",
            fix="Pass --overwrite to overwrite, or choose a different -o prefix",
        )
    if existing and overwrite:
        log.warning("--overwrite: removing %d existing output files", len(existing))
        for p in existing:
            p.unlink()


def _run_stages(
    *,
    args: ExtractCliArgs,
    aadr_df: Any,
    aadr_build: str,
    chain: Path,
    ref_fasta: Path,
    sample_name: str,
    pop_name: str,
    no_lift: bool,
    td_lift: Path,
    td_transform: Path,
    td_call: Path,
) -> tuple[ExtractCounters, rejoin.RejoinOutput]:
    """Execute Stages 1-4 (or just Stage 3 for the no-lift fast path).

    Returns:
        (counters, rejoin_out) — counters is the ExtractCounters aggregate;
        rejoin_out carries the populated sidecar dict + per-variant rows for
        the writers in `_write_outputs`. Returning both lets the orchestrator
        avoid threading the sidecar through ExtractCounters (which would
        couple the JSON-serializable schema to per-sample diagnostics).
    """
    if no_lift:
        snp_path = td_transform / "aadr_native.snp"
        bed_path = td_transform / "aadr_native.bed"
        _write_aadr_native_snp_and_bed(aadr_df, snp_path, bed_path)
        s3_counters = pileup_call.run_pileup_call(
            bam_path=args.bam, snp_path=snp_path, bed_path=bed_path,
            target_fasta_path=ref_fasta,
            output_prefix=td_call / "user_native",
            sample_name=sample_name, pop_name=pop_name,
            seed=args.seed, threads=args.threads,
            min_mapq=args.min_mapq, min_baseq=args.min_baseq,
            no_baq=args.no_baq, no_thread_cap=args.no_thread_cap,
        )
        rejoin_out = rejoin._no_lift_fast_path_finalize(
            pileupcaller_eig_prefix=td_call / "user_native",
            aadr_df=aadr_df,
            output_prefix=args.output_prefix,
            sample_name=sample_name, pop_name=pop_name, sex=args.sex,
            emit_per_variant_rows=(args.report_tsv is not None),
        )
        counters = ExtractCounters(
            stage_1_lift=None, stage_2_transform=None,
            stage_3_call=s3_counters,
            # No Stage 4 for fast path; finalize wallclock folds into total.
            stage_4_rejoin=None,
            coverage=rejoin_out.coverage_counters,
            gates={"liftover_yield": "N/A", "coverage": "PASS"},
            wallclock_total_seconds=0.0,
        )
        return counters, rejoin_out

    sites_vcf_path = td_lift / "aadr_sites.vcf"
    s1_input_filters = sites_vcf.build_sites_vcf(
        aadr_df=aadr_df, output_path=sites_vcf_path,
        aadr_build=aadr_build,  # type: ignore[arg-type]
        palindrome_filter=not args.keep_palindromes,
    )
    lifted_vcf = td_lift / "aadr_lifted.vcf"
    rejected_vcf = td_lift / "aadr_rejected.vcf"
    s1 = lift.lift_aadr_sites(
        sites_vcf_path=sites_vcf_path,
        chain_path=chain, target_fasta_path=ref_fasta,
        output_lifted_vcf=lifted_vcf, output_rejected_vcf=rejected_vcf,
        input_filter_counters=s1_input_filters,
        picard_mem=args.picard_mem, picard_max_records=args.picard_max_records,
        yield_fail_pct=args.liftover_yield_fail_pct,
        yield_warn_pct=args.liftover_yield_warn_pct,
    )

    snp_path = td_transform / "aadr_hg38.snp"
    bed_path = td_transform / "aadr_hg38.bed"
    s2 = transform.build_pileupcaller_snp_and_bed(
        lifted_vcf_path=lifted_vcf,
        output_snp_path=snp_path, output_bed_path=bed_path,
        alt_contig_filter=not args.keep_alt_contigs,
    )

    s3 = pileup_call.run_pileup_call(
        bam_path=args.bam, snp_path=snp_path, bed_path=bed_path,
        target_fasta_path=ref_fasta,
        output_prefix=td_call / "user_hg38",
        sample_name=sample_name, pop_name=pop_name,
        seed=args.seed, threads=args.threads,
        min_mapq=args.min_mapq, min_baseq=args.min_baseq,
        no_baq=args.no_baq, no_thread_cap=args.no_thread_cap,
    )

    rejoin_out = rejoin.rejoin_aadr_frame(
        pileupcaller_eig_prefix=td_call / "user_hg38",
        aadr_df=aadr_df, lifted_vcf_path=lifted_vcf,
        output_prefix=args.output_prefix,
        sample_name=sample_name, pop_name=pop_name, sex=args.sex,
        emit_per_variant_rows=(args.report_tsv is not None),
    )

    counters = ExtractCounters(
        stage_1_lift=s1, stage_2_transform=s2,
        stage_3_call=s3, stage_4_rejoin=rejoin_out.stage_4_counters,
        coverage=rejoin_out.coverage_counters,
        gates={
            "liftover_yield": "WARN" if s1.liftover_yield_warning else "PASS",
            "coverage": "PASS",  # set fully by gate eval
        },
        wallclock_total_seconds=0.0,
    )
    return counters, rejoin_out


def _evaluate_coverage_gate(counters: ExtractCounters, args: ExtractCliArgs) -> None:
    """Apply min/warn coverage thresholds. Raises CoverageGateFailure on min violation."""
    cov = counters.coverage
    n = cov.non_missing_autosomal_calls
    if n < args.min_coverage:
        raise CoverageGateFailure(
            what=f"non-missing autosomal calls {n:,}",
            why=f"below --min-coverage {args.min_coverage:,}",
            fix=(
                "BAM may be low-coverage WGS, exome capture, or wrong panel. "
                "Run `pileup-aadr coverage <bam>` to inspect coverage."
            ),
        )
    if n < args.warn_coverage:
        log.warning(
            "Coverage %d below --warn-coverage %d (proceeding with low-power flag)",
            n, args.warn_coverage,
        )
        counters.coverage = dataclasses.replace(cov, coverage_warning=True)
        counters.gates["coverage"] = "WARN"


def _finalize_counters(
    counters: ExtractCounters, total_wallclock: float
) -> ExtractCounters:
    """Set total wallclock; return updated counters."""
    return dataclasses.replace(counters, wallclock_total_seconds=round(total_wallclock, 1))


def _write_outputs(
    *,
    args: ExtractCliArgs,
    counters: ExtractCounters,
    rejoin_out: rejoin.RejoinOutput,
    bam_format: str,
    bam_build: str,
    aadr_total: int,
    aadr_build: str,
    ref_fasta: Path,
    chain: Path,
    tool_versions: dict[str, str],
) -> None:
    """Sidecar JSON, JSON report, per-variant TSV, stdout summary."""
    sidecar_path = Path(f"{args.output_prefix}.pseudohaploid.json")
    output.write_pseudohaploid_sidecar(sidecar_path, rejoin_out.pseudohaploid_sidecar)

    if args.report_json:
        config_dict = {
            f.name: getattr(args, f.name)
            for f in dataclasses.fields(args)
        }
        config_dict = {
            k: (str(v) if isinstance(v, Path) else v)
            for k, v in config_dict.items()
        }
        input_meta = {
            "bam_path": str(args.bam),
            "bam_format": bam_format,
            "bam_build": bam_build,
            "bam_sample_name_source": _classify_sample_name_source(
                args.sample_name, args.bam
            ),
            "aadr_snp_path": str(args.aadr_snp),
            "aadr_build": aadr_build,
            "aadr_input_rows": aadr_total,
            "chain_path": str(chain),
            "chain_sha256": _read_chain_sha256_for_report(chain),
            "target_fasta_path": str(ref_fasta),
        }
        output_bytes = _compute_output_bytes(args.output_prefix)
        output_meta = {
            "prefix": str(args.output_prefix),
            "geno_bytes": output_bytes["geno"],
            "snp_bytes": output_bytes["snp"],
            "ind_bytes": output_bytes["ind"],
            "pseudohaploid_sidecar": str(sidecar_path),
        }
        output.write_json_report(
            args.report_json, counters,
            config=config_dict, tool_versions=tool_versions,
            input_meta=input_meta, output_meta=output_meta,
        )

    if args.report_tsv:
        output.write_per_variant_tsv(
            args.report_tsv, iter(rejoin_out.per_variant_rows)
        )

    output_bytes = _compute_output_bytes(args.output_prefix)
    output.write_stdout_summary(
        counters, bam_path=args.bam, bam_format=bam_format, bam_build=bam_build,
        bam_coverage=None,
        aadr_path=args.aadr_snp, aadr_total=aadr_total,
        ref_fasta=ref_fasta, chain_path=chain,
        output_prefix=args.output_prefix, output_bytes=output_bytes, sex=args.sex,
    )


def _classify_sample_name_source(explicit_cli: str | None, bam: Path) -> str:
    """Identify which resolution path produced the IID, for JSON report provenance.

    Returns:
        "cli" if --sample-name was explicitly passed; "rg_sm" if BAM @RG SM:
        was the source; "filename" if the IID came from the BAM filename stem.
    """
    if explicit_cli is not None:
        return "cli"
    from .format_detect import _extract_rg_sms
    if _extract_rg_sms(bam):
        return "rg_sm"
    return "filename"


def _read_chain_sha256_for_report(chain_path: Path) -> str:
    """SHA-256 of the chain file used (for reproducibility audit)."""
    return hashlib.sha256(chain_path.read_bytes()).hexdigest()


def _compute_output_bytes(prefix: Path) -> dict[str, int]:
    """Stat the four output artifacts and return byte counts."""
    return {
        "geno": Path(f"{prefix}.geno").stat().st_size,
        "snp": Path(f"{prefix}.snp").stat().st_size,
        "ind": Path(f"{prefix}.ind").stat().st_size,
        "pseudohaploid_json": Path(f"{prefix}.pseudohaploid.json").stat().st_size,
    }


def _write_aadr_native_snp_and_bed(
    aadr_df: Any, snp_path: Path, bed_path: Path,
) -> None:
    """No-lift fast path: write pileupCaller .snp + BED directly from AADR DF.

    Both files use the same encoding as Stage 2's transform.py output (numeric
    chrom in `.snp`, chr-prefixed in BED) for consistency with the lift-path.
    """
    from .format_detect import normalize_chrom
    from .transform import _CANONICAL_CHROM_RE, _CHROM_TO_NUMERIC

    snp_path.parent.mkdir(parents=True, exist_ok=True)
    bed_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with open(snp_path, "w") as snp_out, open(bed_path, "w") as bed_out:
        for rsid, row in aadr_df.iterrows():
            chrom_raw = row["chrom_int"]
            chrom_chr = normalize_chrom(chrom_raw)
            if chrom_chr is None or _CANONICAL_CHROM_RE.match(chrom_chr) is None:
                continue
            chrom_numeric = _CHROM_TO_NUMERIC.get(chrom_chr)
            if chrom_numeric is None:
                continue
            pos = int(row["pos_bp"])
            snp_out.write(
                f"{rsid}\t{chrom_numeric}\t{row['gen_morgans']}\t{pos}\t"
                f"{row['ref']}\t{row['alt']}\n"
            )
            bed_out.write(f"{chrom_chr}\t{pos - 1}\t{pos}\n")
            n_written += 1
    log.info("Wrote no-lift fast-path .snp + BED: %d sites", n_written)


__all__ = ["run_extract"]
