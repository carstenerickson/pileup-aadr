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
    # AADR through v66 is hg19-native, so v0.1 always sees aadr_build="hg19" in practice.
    # The hg38 branch is here for v0.2 forward-compat (no-cost switch table).

    palindrome_drops = 0
    non_snp_drops = 0
    chroms_present: set[str] = set()
    rows_to_write: list[tuple[str, int, str, str, str]] = []

    # Pass 1: filter + collect chrom set
    for rsid, row in aadr_df.iterrows():
        chrom = normalize_chrom(row["chrom_int"])
        if chrom is None or chrom not in chrom_lengths:
            log.debug("Skipping rsid=%s with unrecognized chrom %r", rsid, row["chrom_int"])
            continue
        ref, alt = row["ref"], row["alt"]
        if non_snp_filter and (
            len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT"
        ):
            non_snp_drops += 1
            continue
        if palindrome_filter and (ref, alt) in _PALINDROMES:
            palindrome_drops += 1
            continue
        rows_to_write.append((chrom, int(row["pos_bp"]), str(rsid), ref, alt))
        chroms_present.add(chrom)

    # Sort: chrom (per CHROM_ORDER) then pos
    chrom_index = {c: i for i, c in enumerate(CHROM_ORDER)}
    rows_to_write.sort(key=lambda r: (chrom_index[r[0]], r[1]))

    # Pass 2: write VCF v4.2
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
        for chrom, pos, rsid, ref, alt in rows_to_write:
            out.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\tAADR_RS={rsid}\n")
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
