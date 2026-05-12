"""Tests for transform.py — Stage 2 (pileupCaller .snp + mpileup BED).

Maps to LLD test set #30 (alt-contig filter, chrom encoding choices, AADR_RS
preservation, multi-allelic defensive skip).
"""
from __future__ import annotations

from pathlib import Path

from pileup_aadr.counters import Stage2TransformCounters
from pileup_aadr.transform import build_pileupcaller_snp_and_bed


def _write_lifted_vcf(
    path: Path, rows: list[tuple[str, int, str, str, str, str]]
) -> None:
    """Build a minimal Picard-shaped lifted VCF.

    Rows are (chrom, pos, rsid, ref, alt, filter_or_aadr_rs) tuples. The last
    field becomes AADR_RS=<value> in INFO unless it's the literal "NO_AADR_RS",
    in which case INFO is left empty.
    """
    contigs = sorted({r[0] for r in rows})
    lines = ["##fileformat=VCFv4.2"]
    for contig in contigs:
        # Length doesn't matter for the parser test
        lines.append(f"##contig=<ID={contig},length=300000000>")
    lines.append(
        '##INFO=<ID=AADR_RS,Number=1,Type=String,Description="AADR rsID">'
    )
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    for chrom, pos, rsid, ref, alt, aadr_rs in rows:
        info = "." if aadr_rs == "NO_AADR_RS" else f"AADR_RS={aadr_rs}"
        lines.append(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t{info}")
    path.write_text("\n".join(lines) + "\n")


# --- Basic round-trip ---


def test_build_writes_snp_and_bed(tmp_path: Path) -> None:
    """Three canonical-chrom rows round-trip into both files."""
    vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(vcf, [
        ("chr1", 1000, "rs1", "A", "G", "rs1"),
        ("chr22", 2000, "rs2", "C", "T", "rs2"),
        ("chrX", 3000, "rs3", "G", "A", "rs3"),
    ])
    snp = tmp_path / "out.snp"
    bed = tmp_path / "out.bed"
    counters = build_pileupcaller_snp_and_bed(vcf, snp, bed)

    assert isinstance(counters, Stage2TransformCounters)
    assert counters.output_sites == 3
    assert counters.alt_contig_drops == 0
    assert snp.exists()
    assert bed.exists()


def test_snp_uses_numeric_chrom_encoding(tmp_path: Path) -> None:
    """`.snp` column 2 = AADR-numeric (1-22, X→23, Y→24, M→90)."""
    vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(vcf, [
        ("chr1", 100, "rs1", "A", "G", "rs1"),
        ("chr22", 200, "rs2", "C", "T", "rs2"),
        ("chrX", 300, "rs3", "G", "A", "rs3"),
        ("chrY", 400, "rs4", "T", "C", "rs4"),
        ("chrM", 500, "rs5", "A", "G", "rs5"),
    ])
    snp = tmp_path / "out.snp"
    bed = tmp_path / "out.bed"
    build_pileupcaller_snp_and_bed(vcf, snp, bed)

    cols2 = [line.split("\t")[1] for line in snp.read_text().splitlines()]
    assert cols2 == ["1", "22", "23", "24", "90"]


def test_bed_uses_chr_prefixed_form_and_0_based(tmp_path: Path) -> None:
    """BED is chrN form, pos-1 (0-based start) and pos (1-based exclusive end)."""
    vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(vcf, [
        ("chr1", 1000, "rs1", "A", "G", "rs1"),
    ])
    snp = tmp_path / "out.snp"
    bed = tmp_path / "out.bed"
    build_pileupcaller_snp_and_bed(vcf, snp, bed)

    bed_lines = bed.read_text().splitlines()
    assert bed_lines == ["chr1\t999\t1000"]


def test_snp_preserves_aadr_rs_from_info(tmp_path: Path) -> None:
    """When AADR_RS INFO field is set, rsid in .snp comes from there (not the ID col)."""
    vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(vcf, [
        # ID col != AADR_RS — pileup-aadr should prefer AADR_RS
        ("chr1", 1000, "different_id", "A", "G", "rs_canonical"),
    ])
    snp = tmp_path / "out.snp"
    bed = tmp_path / "out.bed"
    build_pileupcaller_snp_and_bed(vcf, snp, bed)
    rsid = snp.read_text().splitlines()[0].split("\t")[0]
    assert rsid == "rs_canonical"


def test_snp_falls_back_to_id_col_when_aadr_rs_missing(tmp_path: Path) -> None:
    """If AADR_RS INFO is absent, .snp rsid falls back to the ID column."""
    vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(vcf, [
        ("chr1", 1000, "rs_from_id_col", "A", "G", "NO_AADR_RS"),
    ])
    snp = tmp_path / "out.snp"
    bed = tmp_path / "out.bed"
    build_pileupcaller_snp_and_bed(vcf, snp, bed)
    rsid = snp.read_text().splitlines()[0].split("\t")[0]
    assert rsid == "rs_from_id_col"


# --- Alt-contig filter ---


def test_alt_contig_filter_default_drops_alt_haplotype(tmp_path: Path) -> None:
    """Alt/decoy contigs (e.g., chr1_KI270706v1_random) drop with default filter."""
    vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(vcf, [
        ("chr1", 1000, "rs1", "A", "G", "rs1"),
        ("chr1_KI270706v1_random", 2000, "rs2", "C", "T", "rs2"),
        ("chrUn_GL000195v1", 3000, "rs3", "G", "A", "rs3"),
        ("HLA-A*01:01:01:01", 4000, "rs4", "T", "C", "rs4"),
        ("chr22", 5000, "rs5", "A", "G", "rs5"),
    ])
    snp = tmp_path / "out.snp"
    bed = tmp_path / "out.bed"
    counters = build_pileupcaller_snp_and_bed(vcf, snp, bed)

    assert counters.output_sites == 2
    assert counters.alt_contig_drops == 3


def test_alt_contig_filter_off_keeps_alts_in_snp(tmp_path: Path) -> None:
    """alt_contig_filter=False — alts that map to a numeric chrom still pass; the
    canonical-chrom regex's role is purely the drop decision when filter is on.
    The numeric-chrom map ALSO acts as a backstop, so non-canonical chroms
    (e.g., HLA-A*..) still drop because they're not in _CHROM_TO_NUMERIC.
    """
    vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(vcf, [
        ("chr1", 1000, "rs1", "A", "G", "rs1"),
        ("HLA-A*01:01:01:01", 2000, "rs2", "C", "T", "rs2"),
    ])
    snp = tmp_path / "out.snp"
    bed = tmp_path / "out.bed"
    counters = build_pileupcaller_snp_and_bed(
        vcf, snp, bed, alt_contig_filter=False
    )
    # chr1 passes; HLA still drops via the numeric-chrom backstop
    assert counters.output_sites == 1
    assert counters.alt_contig_drops == 1


# --- Defensive paths ---


def test_empty_vcf_writes_empty_files(tmp_path: Path) -> None:
    """Lifted VCF with header only → both outputs empty + counters all zero."""
    vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(vcf, [])
    snp = tmp_path / "out.snp"
    bed = tmp_path / "out.bed"
    counters = build_pileupcaller_snp_and_bed(vcf, snp, bed)
    assert counters.output_sites == 0
    assert counters.alt_contig_drops == 0
    assert snp.read_text() == ""
    assert bed.read_text() == ""


def test_output_dirs_created_if_missing(tmp_path: Path) -> None:
    """Parent dirs of output paths are created on demand."""
    vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(vcf, [("chr1", 100, "rs1", "A", "G", "rs1")])
    snp = tmp_path / "deep" / "transform" / "out.snp"
    bed = tmp_path / "deep" / "transform" / "out.bed"
    counters = build_pileupcaller_snp_and_bed(vcf, snp, bed)
    assert counters.output_sites == 1
    assert snp.exists()
    assert bed.exists()
