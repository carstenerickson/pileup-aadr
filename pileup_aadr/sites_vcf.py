"""Pre-Stage 1: construct a minimal VCF v4.2 from the AADR DataFrame.

Picard's `LiftoverVcf` expects a properly-formed VCF as input — not an EIGENSOFT
`.snp`. This module emits one in the AADR build's coordinate system (`hg19` for
v0.1) so Picard can lift it to the target build.

Why custom Python instead of pysam.VariantFile: pysam's writer goes through htslib's
BCF binary serialization then back to VCF text, ~3x slower than direct text write
at the 1.2M-record scale (verified at design time). The constructed VCF is trivial
enough (no FORMAT, no GT, single INFO field) that the custom writer is simpler and
faster. We use pysam for read-side operations elsewhere.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Final, Literal

import pandas as pd

from .chrom_lengths import CHROM_ORDER, HG19_CHROM_LENGTHS, HG38_CHROM_LENGTHS
from .counters import Stage1InputFilters
from .format_detect import normalize_chrom

log = logging.getLogger(__name__)

# Strand-ambiguous (palindromic) allele pairs: a SNP whose two alleles complement
# each other has identical-looking reads on forward and reverse strands. ~8% of
# biallelic SNPs in 1240k. Picard's RECOVER_SWAPPED_REF_ALT can't tell whether such
# a site is a genuine strand swap or a REF/ALT swap, so we drop them at Stage 1.
_PALINDROMES: Final[frozenset[tuple[str, str]]] = frozenset(
    {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
)
_ACGT: Final[frozenset[str]] = frozenset("ACGT")


def build_sites_vcf(
    aadr_df: pd.DataFrame,
    output_path: Path,
    aadr_build: Literal["hg19", "hg38"],
    palindrome_filter: bool = True,
    non_snp_filter: bool = True,
) -> Stage1InputFilters:
    """Write a minimal VCF v4.2 from the AADR DataFrame.

    Single-base records by construction (palindrome + non-SNP filters drop bad rows
    pre-write). Picard's chain-boundary straddling failure mode (which kills gVCF lifts)
    is structurally impossible here.

    Args:
        aadr_df: from `format_detect.parse_aadr_snp`; index=rsid,
            cols=[chrom_int, gen_morgans, pos_bp, ref, alt]
        output_path: where to write the VCF (typically `<tempdir>/aadr_sites.vcf`)
        aadr_build: "hg19" or "hg38" — the build of AADR's coordinates. The constructed
            VCF describes the SOURCE coordinate system; Picard reads the chain to
            produce TARGET-build output. ##contig declarations match the AADR build.
        palindrome_filter: drop A/T + C/G ambiguous SNPs (default True; ~8% of biallelic)
        non_snp_filter: drop rows where REF or ALT > 1 char or non-ACGT (default True;
            defensive — AADR-shipped files don't have these)

    Returns:
        Stage1InputFilters: counts dropped per filter + rows written.
    """
    chrom_lengths = HG19_CHROM_LENGTHS if aadr_build == "hg19" else HG38_CHROM_LENGTHS

    # Vectorized filter pass — avoids 1.2M pandas.Series allocations from iterrows
    chrom_chr_series = aadr_df["chrom_int"].map(normalize_chrom)
    valid_chrom_mask = chrom_chr_series.notna() & chrom_chr_series.isin(chrom_lengths)

    snp_ok = pd.Series(True, index=aadr_df.index)
    non_snp_drops = 0
    if non_snp_filter:
        snp_ok = (
            (aadr_df["ref"].str.len() == 1) & aadr_df["ref"].isin(_ACGT) &
            (aadr_df["alt"].str.len() == 1) & aadr_df["alt"].isin(_ACGT)
        )
        non_snp_drops = int((valid_chrom_mask & ~snp_ok).sum())

    pal_mask = pd.Series(False, index=aadr_df.index)
    palindrome_drops = 0
    if palindrome_filter:
        pal_mask = (
            ((aadr_df["ref"] == "A") & (aadr_df["alt"] == "T")) |
            ((aadr_df["ref"] == "T") & (aadr_df["alt"] == "A")) |
            ((aadr_df["ref"] == "C") & (aadr_df["alt"] == "G")) |
            ((aadr_df["ref"] == "G") & (aadr_df["alt"] == "C"))
        )
        palindrome_drops = int((valid_chrom_mask & snp_ok & pal_mask).sum())

    write_mask = valid_chrom_mask & snp_ok & ~pal_mask

    # parse_aadr_snp guarantees CHROM_ORDER x pos_bp sort; boolean indexing preserves row order.
    filtered = aadr_df[write_mask].copy()
    filtered["chrom_chr_helper"] = chrom_chr_series[write_mask].values

    chroms_present: set[str] = set(filtered["chrom_chr_helper"].unique())

    rows_written = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        for chrom in CHROM_ORDER:
            if chrom in chroms_present:
                out.write(f"##contig=<ID={chrom},length={chrom_lengths[chrom]}>\n")
        out.write(
            '##INFO=<ID=AADR_RS,Number=1,Type=String,'
            'Description="AADR original SNP name (rsID); preserved through Picard lift">\n'
        )
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in filtered.itertuples():
            chrom = row.chrom_chr_helper
            rsid = row.Index
            out.write(
                f"{chrom}\t{row.pos_bp}\t{rsid}\t{row.ref}\t{row.alt}\t.\tPASS\tAADR_RS={rsid}\n"
            )
            rows_written += 1

    log.info(
        "Built sites VCF at %s: %d rows written, %d palindromes dropped, %d non-SNPs dropped",
        output_path,
        rows_written,
        palindrome_drops,
        non_snp_drops,
    )
    return Stage1InputFilters(
        palindrome_drops=palindrome_drops,
        non_snp_drops=non_snp_drops,
        non_autosome_drops=0,
        rows_written=rows_written,
    )


__all__ = ["build_sites_vcf"]
