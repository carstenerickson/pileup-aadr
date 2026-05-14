"""Tests for rejoin.py — Stage 4 + no-lift fast-path finalizer.

Maps to LLD test set #31 (SwappedAlleles dosage inversion correctness, allele-
mismatch defensive drops, coverage gate fraction = autosomal_calls / autosomal_aadr,
no-lift fast path matches lift path on the equivalence case).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pileup_aadr.errors import PileupAadrInternalError
from pileup_aadr.rejoin import (
    AADR_LOOKUP_FIELDS,
    AadrLookup,
    GENO_HET,
    GENO_HOM_ALT,
    GENO_HOM_REF,
    GENO_MISSING,
    SWAP_DOSAGE,
    RejoinOutput,
    _no_lift_fast_path_finalize,
    build_merged_lookup,
    build_swap_lookup,
    rejoin_aadr_frame,
)

# --- Fixture builders ---


def _aadr_df(rows: list[tuple[str, str, float, int, str, str]]) -> pd.DataFrame:
    """Build an AADR DataFrame indexed by rsid (matches format_detect.parse_aadr_snp output)."""
    df = pd.DataFrame(
        rows,
        columns=["rsid", "chrom_int", "gen_morgans", "pos_bp", "ref", "alt"],
    )
    return df.set_index("rsid", verify_integrity=True)


def _write_pc_triplet(
    prefix: Path,
    rows: list[tuple[str, str, float, int, str, str, str]],
) -> None:
    """Write a synthetic pileupCaller triplet at <prefix>.{geno,snp}.

    Each row: (rsid, chrom_int, gen, pos, ref, alt, geno_char).
    .ind is not written here (rejoin doesn't read it).
    """
    geno_path = Path(f"{prefix}.geno")
    snp_path = Path(f"{prefix}.snp")
    geno_path.parent.mkdir(parents=True, exist_ok=True)
    with open(geno_path, "w") as gh, open(snp_path, "w") as sh:
        for rsid, chrom, gen, pos, ref, alt, geno in rows:
            gh.write(geno + "\n")
            sh.write(f"{rsid}\t{chrom}\t{gen}\t{pos}\t{ref}\t{alt}\n")


def _write_lifted_vcf(
    path: Path, rows: list[tuple[str, int, str, str, str, bool]]
) -> None:
    """Build a Picard-shaped lifted VCF.

    Rows are (chrom, pos, aadr_rs, ref, alt, swapped_alleles) tuples. When
    swapped_alleles=True, the SwappedAlleles flag is set in INFO.
    """
    contigs = sorted({r[0] for r in rows})
    lines = ["##fileformat=VCFv4.2"]
    for contig in contigs:
        lines.append(f"##contig=<ID={contig},length=300000000>")
    lines.append(
        '##INFO=<ID=AADR_RS,Number=1,Type=String,Description="AADR rsID">'
    )
    lines.append(
        '##INFO=<ID=SwappedAlleles,Number=0,Type=Flag,'
        'Description="REF/ALT swapped during lift">'
    )
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    for chrom, pos, aadr_rs, ref, alt, swapped in rows:
        info_parts = [f"AADR_RS={aadr_rs}"]
        if swapped:
            info_parts.append("SwappedAlleles")
        lines.append(
            f"{chrom}\t{pos}\t{aadr_rs}\t{ref}\t{alt}\t.\tPASS\t{';'.join(info_parts)}"
        )
    path.write_text("\n".join(lines) + "\n")


def _make_lookup(aadr: pd.DataFrame, lifted_vcf: Path) -> AadrLookup:
    """Build AadrLookup from a DataFrame + lifted VCF (test helper)."""
    swap = build_swap_lookup(lifted_vcf)
    return build_merged_lookup(aadr, swap)


# --- SWAP_DOSAGE table ---


def test_swap_dosage_table_correctness() -> None:
    """0<->2, 1<->1, 9<->9. Inversion is involutive (swap of swap == identity)."""
    assert SWAP_DOSAGE[GENO_HOM_REF] == GENO_HOM_ALT
    assert SWAP_DOSAGE[GENO_HOM_ALT] == GENO_HOM_REF
    assert SWAP_DOSAGE[GENO_HET] == GENO_HET
    assert SWAP_DOSAGE[GENO_MISSING] == GENO_MISSING
    for k in (GENO_HOM_REF, GENO_HET, GENO_HOM_ALT, GENO_MISSING):
        assert SWAP_DOSAGE[SWAP_DOSAGE[k]] == k


# --- rejoin_aadr_frame: passthrough (no swap) ---


def test_passthrough_no_swap_writes_eigenstrat_triplet(tmp_path: Path) -> None:
    """No SwappedAlleles flag + alleles match → .geno passes through unchanged."""
    aadr = _aadr_df([
        ("rs1", "1", 0.0, 1000, "A", "G"),
        ("rs2", "22", 0.0, 2000, "C", "T"),
    ])
    pc_prefix = tmp_path / "call" / "user_hg38"
    # pileupCaller output lifted to hg38 (different pos but same alleles)
    _write_pc_triplet(pc_prefix, [
        ("rs1", "1", 0.0, 5000, "A", "G", GENO_HOM_REF),
        ("rs2", "22", 0.0, 6000, "C", "T", GENO_HET),
    ])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [
        ("chr1", 5000, "rs1", "A", "G", False),
        ("chr22", 6000, "rs2", "C", "T", False),
    ])

    out_prefix = tmp_path / "out" / "carsten_pseudohaploid"
    result = rejoin_aadr_frame(
        pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "Carsten", "TestPop"
    )

    assert isinstance(result, RejoinOutput)
    assert result.stage_4_counters.output_variants == 2
    assert result.stage_4_counters.ref_alt_swap_count == 0
    assert result.stage_4_counters.allele_mismatch_drops == 0
    geno_lines = (out_prefix.with_suffix(".geno")).read_text().splitlines()
    assert geno_lines == [GENO_HOM_REF, GENO_HET]


def test_output_snp_uses_aadr_hg19_coords(tmp_path: Path) -> None:
    """Output .snp must carry AADR's chrom_int + Morgans + hg19 pos, not hg38."""
    aadr = _aadr_df([
        ("rs1", "1", 0.123, 999_000, "A", "G"),
    ])
    pc_prefix = tmp_path / "call"
    _write_pc_triplet(pc_prefix, [
        ("rs1", "1", 0.0, 5_000_000, "A", "G", GENO_HOM_REF),  # hg38 pos
    ])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [("chr1", 5_000_000, "rs1", "A", "G", False)])

    out_prefix = tmp_path / "out"
    rejoin_aadr_frame(pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "S", "P")

    snp_line = (out_prefix.with_suffix(".snp")).read_text().splitlines()[0]
    parts = snp_line.split("\t")
    assert parts == ["rs1", "1", "0.123", "999000", "A", "G"]


def test_ind_file_written_with_user_sex(tmp_path: Path) -> None:
    """`.ind` is 3-column TSV: IID \\t SEX \\t POP."""
    aadr = _aadr_df([("rs1", "1", 0.0, 1000, "A", "G")])
    pc_prefix = tmp_path / "call"
    _write_pc_triplet(pc_prefix, [("rs1", "1", 0.0, 5000, "A", "G", GENO_HOM_REF)])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [("chr1", 5000, "rs1", "A", "G", False)])

    out_prefix = tmp_path / "out"
    rejoin_aadr_frame(pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "Carsten", "TestPop", sex="M")

    ind = out_prefix.with_suffix(".ind").read_text().rstrip("\n")
    assert ind == "Carsten\tM\tTestPop"


# --- rejoin_aadr_frame: SwappedAlleles dosage inversion ---


def test_swap_inverts_dosage(tmp_path: Path) -> None:
    """SwappedAlleles flag → 0->2, 2->0, 1 stays, 9 stays."""
    # AADR has REF=A, ALT=G — Picard swapped these to REF=G, ALT=A in lifted VCF
    aadr = _aadr_df([
        ("rs_homref", "1", 0.0, 1000, "A", "G"),
        ("rs_het", "1", 0.0, 2000, "A", "G"),
        ("rs_homalt", "1", 0.0, 3000, "A", "G"),
        ("rs_miss", "1", 0.0, 4000, "A", "G"),
    ])
    pc_prefix = tmp_path / "call"
    # pileupCaller called against the SWAPPED reference (REF=G, ALT=A in hg38)
    _write_pc_triplet(pc_prefix, [
        ("rs_homref", "1", 0.0, 5000, "G", "A", GENO_HOM_REF),
        ("rs_het", "1", 0.0, 6000, "G", "A", GENO_HET),
        ("rs_homalt", "1", 0.0, 7000, "G", "A", GENO_HOM_ALT),
        ("rs_miss", "1", 0.0, 8000, "G", "A", GENO_MISSING),
    ])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [
        ("chr1", 5000, "rs_homref", "G", "A", True),
        ("chr1", 6000, "rs_het", "G", "A", True),
        ("chr1", 7000, "rs_homalt", "G", "A", True),
        ("chr1", 8000, "rs_miss", "G", "A", True),
    ])

    out_prefix = tmp_path / "out"
    result = rejoin_aadr_frame(pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "S", "P")

    geno_lines = out_prefix.with_suffix(".geno").read_text().splitlines()
    # 0->2 (swap), 1->1, 2->0 (swap), 9->9
    assert geno_lines == [GENO_HOM_ALT, GENO_HET, GENO_HOM_REF, GENO_MISSING]
    assert result.stage_4_counters.ref_alt_swap_count == 4


# --- Defensive paths ---


def test_swap_flag_but_alleles_dont_swap_drops_with_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """SwappedAlleles flag + alleles don't actually match swap pattern → drop."""
    import logging
    aadr = _aadr_df([("rs1", "1", 0.0, 1000, "A", "G")])
    pc_prefix = tmp_path / "call"
    # Lifted ref/alt don't match either AADR's frame OR a swap
    _write_pc_triplet(pc_prefix, [("rs1", "1", 0.0, 5000, "C", "T", GENO_HOM_REF)])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [("chr1", 5000, "rs1", "C", "T", True)])
    caplog.set_level(logging.WARNING, logger="pileup_aadr.rejoin")

    out_prefix = tmp_path / "out"
    result = rejoin_aadr_frame(pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "S", "P")

    assert result.stage_4_counters.allele_mismatch_drops == 1
    assert result.stage_4_counters.output_variants == 0
    assert any("SwappedAlleles" in r.message for r in caplog.records)


def test_no_swap_flag_but_alleles_differ_drops_with_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """No swap flag but lifted alleles differ from AADR → drop."""
    import logging
    aadr = _aadr_df([("rs1", "1", 0.0, 1000, "A", "G")])
    pc_prefix = tmp_path / "call"
    _write_pc_triplet(pc_prefix, [("rs1", "1", 0.0, 5000, "C", "T", GENO_HOM_REF)])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [("chr1", 5000, "rs1", "C", "T", False)])
    caplog.set_level(logging.WARNING, logger="pileup_aadr.rejoin")

    out_prefix = tmp_path / "out"
    result = rejoin_aadr_frame(pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "S", "P")

    assert result.stage_4_counters.allele_mismatch_drops == 1
    assert result.stage_4_counters.output_variants == 0


def test_rsid_in_pc_output_missing_from_aadr_skipped(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """rsID in pileupCaller output but not in AADR DataFrame → skip with WARNING."""
    import logging
    aadr = _aadr_df([("rs1", "1", 0.0, 1000, "A", "G")])
    pc_prefix = tmp_path / "call"
    _write_pc_triplet(pc_prefix, [
        ("rs1", "1", 0.0, 5000, "A", "G", GENO_HOM_REF),
        ("rs_orphan", "1", 0.0, 6000, "C", "T", GENO_HET),
    ])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [
        ("chr1", 5000, "rs1", "A", "G", False),
        ("chr1", 6000, "rs_orphan", "C", "T", False),
    ])
    caplog.set_level(logging.WARNING, logger="pileup_aadr.rejoin")

    out_prefix = tmp_path / "out"
    result = rejoin_aadr_frame(pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "S", "P")

    assert result.stage_4_counters.rsid_matched == 1
    assert result.stage_4_counters.output_variants == 1
    assert any("missing from AADR" in r.message for r in caplog.records)


def test_malformed_snp_row_raises(tmp_path: Path) -> None:
    """A pileupCaller .snp row with != 6 cols → PileupAadrInternalError."""
    aadr = _aadr_df([("rs1", "1", 0.0, 1000, "A", "G")])
    pc_prefix = tmp_path / "call"
    pc_prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{pc_prefix}.geno").write_text("0\n")
    # Only 4 columns instead of 6
    Path(f"{pc_prefix}.snp").write_text("rs1\t1\t0.0\t5000\n")
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [("chr1", 5000, "rs1", "A", "G", False)])

    out_prefix = tmp_path / "out"
    with pytest.raises(PileupAadrInternalError, match="6 cols"):
        rejoin_aadr_frame(pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "S", "P")


# --- Coverage counters ---


def test_coverage_fraction_uses_autosomes_only(tmp_path: Path) -> None:
    """Coverage = non_missing_autosomal / autosomal_total. Sex chroms in per-chrom but not gated."""
    aadr = _aadr_df([
        ("rs_a1", "1", 0.0, 1000, "A", "G"),
        ("rs_a2", "1", 0.0, 2000, "C", "T"),
        ("rs_a3", "22", 0.0, 3000, "A", "G"),  # autosomal
        ("rs_x", "23", 0.0, 4000, "A", "G"),  # X — not in coverage denominator
    ])
    pc_prefix = tmp_path / "call"
    _write_pc_triplet(pc_prefix, [
        ("rs_a1", "1", 0.0, 5000, "A", "G", GENO_HOM_REF),
        ("rs_a2", "1", 0.0, 6000, "C", "T", GENO_MISSING),  # 1 of 3 autosomes missing
        ("rs_a3", "22", 0.0, 7000, "A", "G", GENO_HET),
        ("rs_x", "23", 0.0, 8000, "A", "G", GENO_HOM_REF),  # X — counted in per-chrom
    ])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [
        ("chr1", 5000, "rs_a1", "A", "G", False),
        ("chr1", 6000, "rs_a2", "C", "T", False),
        ("chr22", 7000, "rs_a3", "A", "G", False),
        ("chrX", 8000, "rs_x", "A", "G", False),
    ])

    out_prefix = tmp_path / "out"
    result = rejoin_aadr_frame(pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "S", "P")

    cov = result.coverage_counters
    assert cov.non_missing_autosomal_calls == 2
    assert cov.coverage_fraction == round(2 / 3, 4)  # 2 calls / 3 autosomal AADR rows
    # X chrom appears in per-chrom (defensive — orchestrator surfaces it)
    assert cov.per_chrom_call_count.get("chrX") == 1


# --- Per-variant rows for --report-tsv ---


def test_emit_per_variant_rows_populates_report(tmp_path: Path) -> None:
    """emit_per_variant_rows=True → per_variant_rows non-empty with action labels."""
    aadr = _aadr_df([
        ("rs_pass", "1", 0.0, 1000, "A", "G"),
        ("rs_swap", "1", 0.0, 2000, "C", "T"),
        ("rs_miss", "1", 0.0, 3000, "A", "G"),
    ])
    pc_prefix = tmp_path / "call"
    _write_pc_triplet(pc_prefix, [
        ("rs_pass", "1", 0.0, 5000, "A", "G", GENO_HOM_REF),
        ("rs_swap", "1", 0.0, 6000, "T", "C", GENO_HOM_REF),  # swapped frame
        ("rs_miss", "1", 0.0, 7000, "A", "G", GENO_MISSING),
    ])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [
        ("chr1", 5000, "rs_pass", "A", "G", False),
        ("chr1", 6000, "rs_swap", "T", "C", True),
        ("chr1", 7000, "rs_miss", "A", "G", False),
    ])

    out_prefix = tmp_path / "out"
    result = rejoin_aadr_frame(
        pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "S", "P",
        emit_per_variant_rows=True,
    )

    actions = [r["action"] for r in result.per_variant_rows]
    assert actions == ["passthrough", "swap", "missing_call"]


# --- Sidecar JSON ---


def test_sidecar_pseudohaploid_classification(tmp_path: Path) -> None:
    """Sidecar marks pseudohaploid=1 with het_count + het_rate populated."""
    aadr = _aadr_df([
        ("rs1", "1", 0.0, 1000, "A", "G"),
        ("rs2", "1", 0.0, 2000, "C", "T"),
    ])
    pc_prefix = tmp_path / "call"
    _write_pc_triplet(pc_prefix, [
        ("rs1", "1", 0.0, 5000, "A", "G", GENO_HET),  # het
        ("rs2", "1", 0.0, 6000, "C", "T", GENO_HOM_REF),
    ])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [
        ("chr1", 5000, "rs1", "A", "G", False),
        ("chr1", 6000, "rs2", "C", "T", False),
    ])

    out_prefix = tmp_path / "out"
    result = rejoin_aadr_frame(pc_prefix, _make_lookup(aadr, lifted_vcf), out_prefix, "Carsten", "P")

    sc = result.pseudohaploid_sidecar
    assert sc["schema_version"] == 1
    sample = sc["samples"]["Carsten"]
    assert sample["pseudohaploid"] == 1
    assert sample["het_count"] == 1
    assert sample["non_missing_autosomal_count"] == 2
    assert sample["het_rate"] == 0.5
    assert sample["calling_mode"] == "randomDiploid"
    assert "no-lift" not in sample["note"]


# --- build_merged_lookup ---


def test_build_merged_lookup_field_layout(tmp_path: Path) -> None:
    """build_merged_lookup returns AadrLookup with correct tuple field order."""
    aadr = _aadr_df([("rs1", "1", 0.123, 999_000, "A", "G")])
    lifted_vcf = tmp_path / "lifted.vcf"
    _write_lifted_vcf(lifted_vcf, [("chr1", 5000, "rs1", "A", "G", False)])
    lookup = _make_lookup(aadr, lifted_vcf)

    assert "rs1" in lookup
    chrom_int, gen_morgans, pos_bp, ref, alt, is_swapped = lookup["rs1"]
    assert chrom_int == "1"
    assert gen_morgans == pytest.approx(0.123)
    assert pos_bp == 999_000
    assert ref == "A"
    assert alt == "G"
    assert is_swapped is False


def test_build_merged_lookup_swap_default_false(tmp_path: Path) -> None:
    """rsID not in lifted VCF → is_swapped defaults to False."""
    aadr = _aadr_df([("rs_absent", "22", 0.0, 50_000, "C", "T")])
    lifted_vcf = tmp_path / "lifted.vcf"
    # Write a VCF with a *different* rsID so rs_absent is absent from swap_lookup
    _write_lifted_vcf(lifted_vcf, [("chr1", 1000, "rs_other", "A", "G", True)])
    lookup = _make_lookup(aadr, lifted_vcf)

    assert "rs_absent" in lookup
    assert lookup["rs_absent"][5] is False  # is_swapped


def test_build_merged_lookup_empty_swap_lookup(tmp_path: Path) -> None:
    """build_merged_lookup with an empty swap_lookup → all is_swapped=False."""
    from pileup_aadr.rejoin import build_merged_lookup
    aadr = _aadr_df([
        ("rs1", "1", 0.0, 1000, "A", "G"),
        ("rs2", "2", 0.0, 2000, "C", "T"),
    ])
    lookup = build_merged_lookup(aadr, {})
    assert all(v[5] is False for v in lookup.values())
    assert len(lookup) == 2


# --- No-lift fast path ---


def test_no_lift_fast_path_copies_triplet(tmp_path: Path) -> None:
    """No-lift path: pileupCaller output IS the AADR-frame triplet, just copied."""
    aadr = _aadr_df([
        ("rs1", "1", 0.0, 1000, "A", "G"),
        ("rs2", "22", 0.0, 2000, "C", "T"),
    ])
    pc_prefix = tmp_path / "call"
    # In the no-lift case, pileupCaller's coords already match AADR
    _write_pc_triplet(pc_prefix, [
        ("rs1", "1", 0.0, 1000, "A", "G", GENO_HOM_REF),
        ("rs2", "22", 0.0, 2000, "C", "T", GENO_HET),
    ])

    out_prefix = tmp_path / "out"
    result = _no_lift_fast_path_finalize(
        pc_prefix, aadr, out_prefix, "Carsten", "TestPop"
    )

    geno_lines = out_prefix.with_suffix(".geno").read_text().splitlines()
    snp_lines = out_prefix.with_suffix(".snp").read_text().splitlines()
    assert geno_lines == [GENO_HOM_REF, GENO_HET]
    assert len(snp_lines) == 2
    assert result.stage_4_counters.output_variants == 2
    assert result.stage_4_counters.ref_alt_swap_count == 0
    assert result.stage_4_counters.allele_mismatch_drops == 0


def test_no_lift_fast_path_sidecar_notes_no_lift(tmp_path: Path) -> None:
    """No-lift sidecar's note mentions the fast path explicitly."""
    aadr = _aadr_df([("rs1", "1", 0.0, 1000, "A", "G")])
    pc_prefix = tmp_path / "call"
    _write_pc_triplet(pc_prefix, [("rs1", "1", 0.0, 1000, "A", "G", GENO_HOM_REF)])

    out_prefix = tmp_path / "out"
    result = _no_lift_fast_path_finalize(pc_prefix, aadr, out_prefix, "S", "P")
    assert "no-lift" in result.pseudohaploid_sidecar["samples"]["S"]["note"]


def test_no_lift_overrides_pc_ind_with_user_sex(tmp_path: Path) -> None:
    """No-lift writes new .ind with user-supplied sex (pileupCaller writes SEX=U)."""
    aadr = _aadr_df([("rs1", "1", 0.0, 1000, "A", "G")])
    pc_prefix = tmp_path / "call"
    _write_pc_triplet(pc_prefix, [("rs1", "1", 0.0, 1000, "A", "G", GENO_HOM_REF)])
    # PileupCaller's .ind would say SEX=U; we should override
    Path(f"{pc_prefix}.ind").write_text("PCSample\tU\tPCPop\n")

    out_prefix = tmp_path / "out"
    _no_lift_fast_path_finalize(
        pc_prefix, aadr, out_prefix, "Carsten", "TestPop", sex="F"
    )
    assert out_prefix.with_suffix(".ind").read_text().rstrip("\n") == "Carsten\tF\tTestPop"
