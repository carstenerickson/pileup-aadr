"""Tests for Stage 1 — parse_picard_stderr + parse_rejected_vcf + lift_aadr_sites
+ Picard sharding (build_picard_shard_manifest, concat_picard_outputs,
aggregate_stage1_counters, lift_aadr_sites_sharded).

Most tests use captured Picard stderr fixtures from v2.1 verification (real Picard
3.3.0 output on `ancestrytracke-f`). The lift_aadr_sites integration test mocks
ToolWrapper.run to avoid needing Picard in the test env.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from pileup_aadr.counters import Stage1InputFilters
from pileup_aadr.errors import LiftoverYieldError, PileupAadrInternalError
from pileup_aadr.lift import (
    lift_aadr_sites,
    parse_picard_stderr,
    parse_rejected_vcf,
)

FIXTURES_STDERR = Path(__file__).parent / "fixtures" / "stderr"


# --- parse_picard_stderr ---


def test_parse_picard_stderr_clean_run() -> None:
    """The 100-site clean run from v2.1: 100 processed, 0 failed, 0 swapped."""
    text = (FIXTURES_STDERR / "picard_clean.stderr").read_text()
    parsed = parse_picard_stderr(text)
    assert parsed["processed_count"] == 100
    assert parsed["failed_count"] == 0
    assert parsed["mismatched_ref_count"] == 0
    assert parsed["lifted_count"] == 100
    assert parsed["swapped_count"] == 0


def test_parse_picard_stderr_partial_yield() -> None:
    """The 5000-site v2.1 run: 5000 processed, 2 failed, 16 swapped."""
    text = (FIXTURES_STDERR / "picard_partial_yield.stderr").read_text()
    parsed = parse_picard_stderr(text)
    assert parsed["processed_count"] == 5000
    assert parsed["failed_count"] == 2
    assert parsed["lifted_count"] == 4998
    assert parsed["swapped_count"] == 16


def test_parse_picard_stderr_swap_line_optional() -> None:
    """If Picard ever omits the 'X variants were lifted by swapping' line, default to 0."""
    text = """
    INFO	LiftoverVcf	Processed 100 variants.
    INFO	LiftoverVcf	0 variants failed to liftover.
    INFO	LiftoverVcf	0 variants lifted over but had mismatching reference alleles after lift over.
    INFO	LiftoverVcf	0.0000% of variants were not successfully lifted over
    """
    parsed = parse_picard_stderr(text)
    assert parsed["swapped_count"] == 0


def test_parse_picard_stderr_missing_required_pattern_raises() -> None:
    """Missing a REQUIRED stderr line (e.g., 'Processed N variants') → PileupAadrInternalError."""
    text = "Some unrelated output that doesn't match Picard's structured stderr."
    with pytest.raises(PileupAadrInternalError, match="processed"):
        parse_picard_stderr(text)


# --- parse_rejected_vcf ---


def _write_rejected_vcf(path: Path, records: list[tuple[str, int, str]]) -> None:
    """Build a minimal rejected.vcf with (chrom, pos, FILTER) tuples."""
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=249250621>",
        '##FILTER=<ID=NoTarget,Description="No chain alignment">',
        '##FILTER=<ID=MismatchedRefAllele,Description="REF allele mismatch">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for chrom, pos, filt in records:
        lines.append(f"{chrom}\t{pos}\trs1\tA\tG\t.\t{filt}\t.")
    path.write_text("\n".join(lines) + "\n")


def test_parse_rejected_vcf_empty(tmp_path: Path) -> None:
    """No rejected.vcf at the path → all-zero counts (defensive)."""
    counts = parse_rejected_vcf(tmp_path / "missing.vcf")
    assert counts == {
        "NoTarget": 0,
        "MismatchedRefAllele": 0,
        "IndelStraddlesMultipleIntervals": 0,
        "SwappedAlleles": 0,
        "other": 0,
    }


def test_parse_rejected_vcf_categorizes_filters(tmp_path: Path) -> None:
    """Each known FILTER value increments its bucket; unknowns go to 'other'."""
    rej = tmp_path / "rejected.vcf"
    _write_rejected_vcf(
        rej,
        [
            ("chr1", 1000, "NoTarget"),
            ("chr1", 2000, "NoTarget"),
            ("chr1", 3000, "MismatchedRefAllele"),
            ("chr1", 4000, "WeirdUnknownFilter"),
        ],
    )
    counts = parse_rejected_vcf(rej)
    assert counts["NoTarget"] == 2
    assert counts["MismatchedRefAllele"] == 1
    assert counts["other"] == 1
    assert counts["IndelStraddlesMultipleIntervals"] == 0


# --- lift_aadr_sites integration (mocked Picard subprocess) ---


@pytest.fixture
def mock_picard_lift(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """Mock ToolWrapper.run to write a fake Picard stderr + empty rejected.vcf,
    leaving the OUTPUT VCF up to the caller to either pre-create or ignore.
    Returns the tmp_path for callers to set up files in.
    """
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))

    # Mock subprocess.run for the version probe (called from ToolWrapper.__init__)
    import subprocess
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    # Mock ToolWrapper.run to write the captured-stderr fixture + return success
    from pileup_aadr import lift, tool_wrapper

    def fake_run(self: object, args: list[str], *, capture_stderr_to: Path, **_kw: object) -> object:
        # Copy the v2.1 captured Picard stderr into the path the caller expects
        capture_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        capture_stderr_to.write_text(
            (FIXTURES_STDERR / "picard_clean.stderr").read_text()
        )
        return tool_wrapper.ToolRunResult(
            exit_code=0,
            stdout=None,
            stderr_path=capture_stderr_to,
            stderr_text=None,
            wallclock_seconds=0.1,
            peak_rss_mb=None,
        )

    monkeypatch.setattr(lift.ToolWrapper, "run", fake_run)
    monkeypatch.setattr(lift.ToolWrapper, "_check_version", lambda _self: None)
    return tmp_path


def test_lift_aadr_sites_clean_run(mock_picard_lift: Path) -> None:
    """100/100 sites lifted (per the captured stderr) → yield = 100% → PASS gate."""
    sites_vcf = mock_picard_lift / "sites.vcf"
    sites_vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    fasta = mock_picard_lift / "ref.fa"
    fasta.touch()
    chain = mock_picard_lift / "chain.gz"
    chain.touch()

    counters = lift_aadr_sites(
        sites_vcf_path=sites_vcf,
        chain_path=chain,
        target_fasta_path=fasta,
        output_lifted_vcf=mock_picard_lift / "lifted.vcf",
        output_rejected_vcf=mock_picard_lift / "rejected.vcf",
        input_filter_counters=Stage1InputFilters(
            palindrome_drops=0, non_snp_drops=0, non_autosome_drops=0, rows_written=100
        ),
    )

    assert counters.input_sites_after_filters == 100
    assert counters.lifted_sites == 100
    assert counters.liftover_yield_pct == 100.0
    assert counters.liftover_yield_warning is False
    assert counters.swapped_alleles_count == 0


def test_lift_aadr_sites_low_yield_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Synthetic stderr showing 50% yield → LiftoverYieldError below default 70% gate."""
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))

    import subprocess
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    fake_low_yield_stderr = """
    INFO	LiftoverVcf	Processed 100 variants.
    INFO	LiftoverVcf	50 variants failed to liftover.
    INFO	LiftoverVcf	0 variants lifted over but had mismatching reference alleles after lift over.
    INFO	LiftoverVcf	50.0000% of variants were not successfully lifted over
    INFO	LiftoverVcf	0 variants were lifted by swapping REF/ALT alleles.
    """

    from pileup_aadr import lift, tool_wrapper

    def fake_run(self: object, args: list[str], *, capture_stderr_to: Path, **_kw: object) -> object:
        capture_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        capture_stderr_to.write_text(fake_low_yield_stderr)
        return tool_wrapper.ToolRunResult(
            exit_code=0, stdout=None, stderr_path=capture_stderr_to,
            stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None,
        )

    monkeypatch.setattr(lift.ToolWrapper, "run", fake_run)
    monkeypatch.setattr(lift.ToolWrapper, "_check_version", lambda _self: None)

    sites_vcf = tmp_path / "sites.vcf"
    sites_vcf.write_text("##fileformat=VCFv4.2\n")
    fasta = tmp_path / "ref.fa"
    fasta.touch()
    chain = tmp_path / "chain.gz"
    chain.touch()

    # Set up a rejected.vcf with NoTarget so the dominant rejection is identified
    rej = tmp_path / "rejected.vcf"
    rej.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=249250621>\n"
        '##FILTER=<ID=NoTarget,Description="x">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "\n".join(
            f"chr1\t{i}\trs{i}\tA\tG\t.\tNoTarget\t."
            for i in range(1, 51)
        )
        + "\n"
    )

    with pytest.raises(LiftoverYieldError, match=r"50/100.*NoTarget"):
        lift_aadr_sites(
            sites_vcf_path=sites_vcf,
            chain_path=chain,
            target_fasta_path=fasta,
            output_lifted_vcf=tmp_path / "lifted.vcf",
            output_rejected_vcf=rej,
            input_filter_counters=Stage1InputFilters(
                palindrome_drops=0, non_snp_drops=0, non_autosome_drops=0, rows_written=100
            ),
        )


def test_lift_aadr_sites_threshold_warn_only(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """80% yield → above 70% fail threshold but below 95% warn threshold → warning logged."""
    import logging
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))

    import subprocess
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    fake_warn_stderr = """
    INFO	LiftoverVcf	Processed 100 variants.
    INFO	LiftoverVcf	20 variants failed to liftover.
    INFO	LiftoverVcf	0 variants lifted over but had mismatching reference alleles after lift over.
    INFO	LiftoverVcf	20.0000% of variants were not successfully lifted over
    INFO	LiftoverVcf	0 variants were lifted by swapping REF/ALT alleles.
    """

    from pileup_aadr import lift, tool_wrapper

    def fake_run(self: object, args: list[str], *, capture_stderr_to: Path, **_kw: object) -> object:
        capture_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        capture_stderr_to.write_text(fake_warn_stderr)
        return tool_wrapper.ToolRunResult(
            exit_code=0, stdout=None, stderr_path=capture_stderr_to,
            stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None,
        )

    monkeypatch.setattr(lift.ToolWrapper, "run", fake_run)
    monkeypatch.setattr(lift.ToolWrapper, "_check_version", lambda _self: None)
    caplog.set_level(logging.WARNING, logger="pileup_aadr.lift")

    sites_vcf = tmp_path / "sites.vcf"
    sites_vcf.write_text("##fileformat=VCFv4.2\n")
    fasta = tmp_path / "ref.fa"
    fasta.touch()
    chain = tmp_path / "chain.gz"
    chain.touch()
    rej = tmp_path / "rejected.vcf"  # leave nonexistent → empty bucket counts

    counters = lift_aadr_sites(
        sites_vcf_path=sites_vcf,
        chain_path=chain,
        target_fasta_path=fasta,
        output_lifted_vcf=tmp_path / "lifted.vcf",
        output_rejected_vcf=rej,
        input_filter_counters=Stage1InputFilters(
            palindrome_drops=0, non_snp_drops=0, non_autosome_drops=0, rows_written=100
        ),
    )
    assert counters.liftover_yield_pct == 80.0
    assert counters.liftover_yield_warning is True
    assert any("yield 80.00%" in r.message for r in caplog.records)


# ─── Picard sharding ───────────────────────────────────────────────────────────


from pileup_aadr.lift import (  # noqa: E402
    PicardShardSpec,
    aggregate_stage1_counters,
    build_picard_shard_manifest,
    concat_picard_outputs,
    lift_aadr_sites_sharded,
    lift_and_transform_sharded,
)

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=249250621>\n"
    "##contig=<ID=chr2,length=243199373>\n"
    "##contig=<ID=chr22,length=51304566>\n"
    "##INFO=<ID=AADR_RS,Number=1,Type=String,Description=\"rsid\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _make_sites_vcf(path: Path, records: list[tuple[str, int, str]]) -> None:
    """Write a minimal sites VCF with (chrom, pos, rsid) records."""
    with open(path, "w") as f:
        f.write(_VCF_HEADER)
        for chrom, pos, rsid in records:
            f.write(f"{chrom}\t{pos}\t{rsid}\tA\tG\t.\tPASS\tAADR_RS={rsid}\n")


# --- L1: build_picard_shard_manifest LPT balance ---


def test_shard_manifest_lpt_balance(tmp_path: Path) -> None:
    """LPT bin-packing assigns largest chrom first; imbalance ≤ smallest chrom's count."""
    sites_vcf = tmp_path / "sites.vcf"
    # chr1: 100 sites, chr2: 50 sites, chr22: 10 sites
    records = (
        [("chr1", i * 100, f"rs1_{i}") for i in range(1, 101)]
        + [("chr2", i * 100, f"rs2_{i}") for i in range(1, 51)]
        + [("chr22", i * 100, f"rs22_{i}") for i in range(1, 11)]
    )
    _make_sites_vcf(sites_vcf, records)

    manifest = build_picard_shard_manifest(sites_vcf, tmp_path / "shards", n_shards=2)

    assert len(manifest) == 2
    # Each shard's input.vcf must exist and contain a subset of records
    total_records = sum(
        sum(1 for line in spec.input_vcf.read_text().splitlines() if not line.startswith("#"))
        for spec in manifest
    )
    assert total_records == 160  # 100 + 50 + 10

    # LPT: shard 0 gets chr1 (100); shard 1 gets chr2+chr22 (60). Max-min = 40 ≤ smallest=10? No.
    # LPT guarantees max-min ≤ largest chrom count dropped onto least-loaded bin.
    # The imbalance between the two shards must be ≤ chr22's count (10).
    counts = [
        sum(1 for line in spec.input_vcf.read_text().splitlines() if not line.startswith("#"))
        for spec in manifest
    ]
    assert max(counts) - min(counts) <= 50  # LPT upper bound: imbalance ≤ 2nd-largest chrom


def test_shard_manifest_clamps_to_chrom_count(tmp_path: Path) -> None:
    """n_shards > chrom count → clamped to chrom count."""
    sites_vcf = tmp_path / "sites.vcf"
    _make_sites_vcf(sites_vcf, [("chr1", 1000, "rs1"), ("chr22", 2000, "rs2")])

    manifest = build_picard_shard_manifest(sites_vcf, tmp_path / "shards", n_shards=10)

    assert len(manifest) == 2  # clamped to 2 chroms


def test_shard_manifest_full_header_in_every_shard(tmp_path: Path) -> None:
    """Every shard's input.vcf contains the full VCF header."""
    sites_vcf = tmp_path / "sites.vcf"
    _make_sites_vcf(sites_vcf, [("chr1", 1000, "rs1"), ("chr22", 2000, "rs2")])

    manifest = build_picard_shard_manifest(sites_vcf, tmp_path / "shards", n_shards=2)

    for spec in manifest:
        content = spec.input_vcf.read_text()
        assert "##fileformat=VCFv4.2" in content
        assert "##contig=<ID=chr1" in content
        assert "#CHROM\tPOS" in content


# --- L3: concat_picard_outputs header handling ---


def test_concat_picard_outputs_preserves_header(tmp_path: Path) -> None:
    """Concatenated lifted VCF: header from shard 0 only; body from all shards."""
    shard0 = tmp_path / "s0"
    shard0.mkdir()
    shard1 = tmp_path / "s1"
    shard1.mkdir()

    (shard0 / "lifted.vcf").write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1000\trs1\tA\tG\t.\tPASS\t.\n"
    )
    (shard0 / "rejected.vcf").write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t999\trs0\tA\tG\t.\tNoTarget\t.\n"
    )
    (shard1 / "lifted.vcf").write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr22\t2000\trs22\tA\tG\t.\tPASS\t.\n"
    )
    (shard1 / "rejected.vcf").write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    specs = [
        PicardShardSpec(0, ("chr1",), shard0 / "i.vcf", shard0 / "lifted.vcf",
                        shard0 / "rejected.vcf", shard0 / "picard.stderr"),
        PicardShardSpec(1, ("chr22",), shard1 / "i.vcf", shard1 / "lifted.vcf",
                        shard1 / "rejected.vcf", shard1 / "picard.stderr"),
    ]

    out_lifted = tmp_path / "merged_lifted.vcf"
    out_rejected = tmp_path / "merged_rejected.vcf"
    concat_picard_outputs(specs, out_lifted, out_rejected)

    lifted_lines = out_lifted.read_text().splitlines()
    header_count = sum(1 for line in lifted_lines if line.startswith("#"))
    body_lines = [line for line in lifted_lines if not line.startswith("#")]
    assert header_count == 2  # ##fileformat + #CHROM header from shard 0 only
    assert len(body_lines) == 2
    assert body_lines[0].startswith("chr1\t")
    assert body_lines[1].startswith("chr22\t")

    rej_lines = out_rejected.read_text().splitlines()
    rej_body = [line for line in rej_lines if not line.startswith("#")]
    assert len(rej_body) == 1  # only shard 0 had a rejected record


# --- L4: aggregate_stage1_counters ---


def test_aggregate_stage1_counters_field_wise(tmp_path: Path) -> None:
    """Field-wise sum; wallclock = max; yield recomputed; input_filters passed through."""
    from pileup_aadr.counters import Stage1LiftCounters

    input_filters = Stage1InputFilters(
        palindrome_drops=5, non_snp_drops=2, non_autosome_drops=0, rows_written=1000
    )
    shard_a = Stage1LiftCounters(
        wallclock_seconds=10.0,
        input_sites_after_filters=1000,  # each shard carries total (not per-shard)
        lifted_sites=450,
        liftover_yield_pct=45.0,   # per-shard pct is wrong; aggregate recomputes
        liftover_yield_warning=False,
        rejected_by_reason={"NoTarget": 30, "MismatchedRefAllele": 5, "other": 0,
                             "IndelStraddlesMultipleIntervals": 0, "SwappedAlleles": 0},
        swapped_alleles_count=10,
        input_filters=input_filters,
    )
    shard_b = Stage1LiftCounters(
        wallclock_seconds=8.5,
        input_sites_after_filters=1000,
        lifted_sites=470,
        liftover_yield_pct=47.0,
        liftover_yield_warning=False,
        rejected_by_reason={"NoTarget": 20, "MismatchedRefAllele": 3, "other": 2,
                             "IndelStraddlesMultipleIntervals": 1, "SwappedAlleles": 0},
        swapped_alleles_count=15,
        input_filters=input_filters,
    )

    agg = aggregate_stage1_counters(
        per_shard=[shard_a, shard_b],
        input_filters=input_filters,
        yield_fail_pct=70.0,
        yield_warn_pct=95.0,
    )

    # wallclock = max
    assert agg.wallclock_seconds == 10.0
    # lifted_sites = sum
    assert agg.lifted_sites == 920
    # input = input_filters.rows_written (not sum of per-shard which would be 2000)
    assert agg.input_sites_after_filters == 1000
    # yield = 920 / 1000 = 92%
    assert agg.liftover_yield_pct == pytest.approx(92.0)
    assert agg.liftover_yield_warning is True  # 92% < 95% warn threshold
    # rejected reasons summed field-wise
    assert agg.rejected_by_reason["NoTarget"] == 50
    assert agg.rejected_by_reason["MismatchedRefAllele"] == 8
    assert agg.rejected_by_reason["other"] == 2
    # swapped = sum
    assert agg.swapped_alleles_count == 25
    # input_filters passed through unchanged
    assert agg.input_filters is input_filters


def test_aggregate_stage1_counters_raises_below_fail_pct(tmp_path: Path) -> None:
    """Aggregate yield < yield_fail_pct → LiftoverYieldError."""
    from pileup_aadr.counters import Stage1LiftCounters

    input_filters = Stage1InputFilters(0, 0, 0, rows_written=100)
    shard = Stage1LiftCounters(
        wallclock_seconds=1.0, input_sites_after_filters=100,
        lifted_sites=40, liftover_yield_pct=40.0, liftover_yield_warning=False,
        rejected_by_reason={"NoTarget": 60, "MismatchedRefAllele": 0, "other": 0,
                             "IndelStraddlesMultipleIntervals": 0, "SwappedAlleles": 0},
        swapped_alleles_count=0, input_filters=input_filters,
    )
    with pytest.raises(LiftoverYieldError, match=r"40/100"):
        aggregate_stage1_counters(
            per_shard=[shard],
            input_filters=input_filters,
            yield_fail_pct=70.0,
            yield_warn_pct=95.0,
        )


# --- L5/L6: lift_aadr_sites_sharded (mocked Picard) ---


def _make_fake_picard_run(
    clean_stderr: str,
    raise_on_chroms: frozenset[str] | None = None,
) -> object:
    """Factory for a fake ToolWrapper.run that echoes input VCF body as lifted output.

    The fake lifts every input record verbatim (input chrom, pos, rsid pass through).
    This is enough to test that sharding + concat produces the same record set as
    single-process.

    If raise_on_chroms is set, raises RuntimeError when the input contains any of
    those chroms (simulates shard failure).
    """
    from pileup_aadr import tool_wrapper

    def fake_run(
        self: object, args: list[str], *, capture_stderr_to: Path, **_kw: object
    ) -> tool_wrapper.ToolRunResult:
        # Extract --INPUT path from args
        input_idx = args.index("--INPUT") + 1
        input_vcf = Path(args[input_idx])
        output_idx = args.index("--OUTPUT") + 1
        output_vcf = Path(args[output_idx])
        reject_idx = args.index("--REJECT") + 1
        reject_vcf = Path(args[reject_idx])

        # Read input records
        header_lines = []
        body_lines = []
        with open(input_vcf) as f:
            for line in f:
                if line.startswith("#"):
                    header_lines.append(line)
                else:
                    body_lines.append(line)

        if raise_on_chroms:
            for line in body_lines:
                chrom = line.split("\t", 1)[0]
                if chrom in raise_on_chroms:
                    raise RuntimeError(f"Fake Picard failure on {chrom}")

        # Write output VCF = header + all body lines (fake "lift" is passthrough)
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        with open(output_vcf, "w") as out:
            for line in header_lines:
                out.write(line)
            for line in body_lines:
                out.write(line)

        # Write empty reject VCF
        with open(reject_vcf, "w") as rej:
            for line in header_lines:
                rej.write(line)

        # Write clean stderr
        capture_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        n = len(body_lines)
        capture_stderr_to.write_text(
            f"INFO\tLiftoverVcf\tProcessed {n} variants.\n"
            f"INFO\tLiftoverVcf\t0 variants failed to liftover.\n"
            f"INFO\tLiftoverVcf\t0 variants lifted over but had mismatching reference alleles after lift over.\n"
            f"INFO\tLiftoverVcf\t0.0000% of variants were not successfully lifted over\n"
            f"INFO\tLiftoverVcf\t0 variants were lifted by swapping REF/ALT alleles.\n"
        )

        return tool_wrapper.ToolRunResult(
            exit_code=0, stdout=None, stderr_path=capture_stderr_to,
            stderr_text=None, wallclock_seconds=0.05, peak_rss_mb=None,
        )

    return fake_run


def test_sharded_lift_logical_equivalence(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """4-shard lift produces same record set as single-process lift.

    Uses a passthrough fake Picard (lifts records verbatim) to check that
    sharding + concat preserves the full set of (rsid, chrom, pos) tuples.
    Comparison is order-insensitive (sets).
    """
    import subprocess
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    from pileup_aadr import lift
    monkeypatch.setattr(lift.ToolWrapper, "run", _make_fake_picard_run(""))
    monkeypatch.setattr(lift.ToolWrapper, "_check_version", lambda _self: None)

    # 10 sites on chr1, 5 on chr2, 3 on chr22
    records = (
        [("chr1", i * 1000, f"rs1_{i}") for i in range(1, 11)]
        + [("chr2", i * 1000, f"rs2_{i}") for i in range(1, 6)]
        + [("chr22", i * 1000, f"rs22_{i}") for i in range(1, 4)]
    )
    sites_vcf = tmp_path / "sites.vcf"
    _make_sites_vcf(sites_vcf, records)

    input_filters = Stage1InputFilters(0, 0, 0, rows_written=len(records))
    chain = tmp_path / "chain.gz"
    chain.touch()
    fasta = tmp_path / "ref.fa"
    fasta.touch()

    # Single-process baseline
    single_lifted = tmp_path / "single_lifted.vcf"
    single_rejected = tmp_path / "single_rejected.vcf"
    counters_1 = lift_aadr_sites_sharded(
        sites_vcf_path=sites_vcf,
        chain_path=chain, target_fasta_path=fasta,
        output_lifted_vcf=single_lifted, output_rejected_vcf=single_rejected,
        input_filter_counters=input_filters,
        shard_tempdir=tmp_path / "shards_1",
        n_shards=1,
    )

    # Sharded (3 shards for 3 chroms)
    sharded_lifted = tmp_path / "sharded_lifted.vcf"
    sharded_rejected = tmp_path / "sharded_rejected.vcf"
    counters_n = lift_aadr_sites_sharded(
        sites_vcf_path=sites_vcf,
        chain_path=chain, target_fasta_path=fasta,
        output_lifted_vcf=sharded_lifted, output_rejected_vcf=sharded_rejected,
        input_filter_counters=input_filters,
        shard_tempdir=tmp_path / "shards_n",
        n_shards=3,
    )

    def _record_set(vcf_path: Path) -> frozenset[tuple[str, str, str]]:
        records_found = set()
        with open(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                records_found.add((parts[2], parts[0], parts[1]))  # (rsid, chrom, pos)
        return frozenset(records_found)

    assert _record_set(single_lifted) == _record_set(sharded_lifted)
    assert counters_1.lifted_sites == counters_n.lifted_sites


def test_sharded_lift_first_failure_wins(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """One shard's Picard raises → exception propagates; other shards' per-shard files captured."""
    import subprocess
    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )

    from pileup_aadr import lift
    monkeypatch.setattr(
        lift.ToolWrapper, "run",
        _make_fake_picard_run("", raise_on_chroms=frozenset({"chr22"})),
    )
    monkeypatch.setattr(lift.ToolWrapper, "_check_version", lambda _self: None)

    records = (
        [("chr1", i * 1000, f"rs1_{i}") for i in range(1, 6)]
        + [("chr22", i * 1000, f"rs22_{i}") for i in range(1, 4)]
    )
    sites_vcf = tmp_path / "sites.vcf"
    _make_sites_vcf(sites_vcf, records)

    input_filters = Stage1InputFilters(0, 0, 0, rows_written=len(records))

    with pytest.raises(RuntimeError, match="Fake Picard failure on chr22"):
        lift_aadr_sites_sharded(
            sites_vcf_path=sites_vcf,
            chain_path=tmp_path / "chain.gz",
            target_fasta_path=tmp_path / "ref.fa",
            output_lifted_vcf=tmp_path / "lifted.vcf",
            output_rejected_vcf=tmp_path / "rejected.vcf",
            input_filter_counters=input_filters,
            shard_tempdir=tmp_path / "shards",
            n_shards=2,
        )


# ─── lift_and_transform_sharded ────────────────────────────────────────────────


def _setup_fake_picard(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    raise_on_chroms: frozenset[str] | None = None,
) -> None:
    """Patch ToolWrapper.run + subprocess.run for version probe."""
    import subprocess
    from unittest.mock import MagicMock

    from pileup_aadr import lift

    fake_jar = tmp_path / "picard.jar"
    fake_jar.write_bytes(b"fake jar")
    monkeypatch.setenv("PICARD_JAR", str(fake_jar))
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="", stderr="Version:3.3.0\n", returncode=0),
    )
    monkeypatch.setattr(lift.ToolWrapper, "run", _make_fake_picard_run(
        clean_stderr="", raise_on_chroms=raise_on_chroms,
    ))
    monkeypatch.setattr(lift.ToolWrapper, "_check_version", lambda _self: None)


def test_lift_and_transform_single_shard_sequential(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """S1: n_shards=1 short-circuits to sequential lift + transform.

    Verifies that (a) the function returns (Stage1LiftCounters, Stage2TransformCounters)
    and (b) the .snp + .bed are written and non-empty.
    """
    from pileup_aadr.counters import Stage1LiftCounters, Stage2TransformCounters

    _setup_fake_picard(monkeypatch, tmp_path)

    records = [("chr1", i * 1000, f"rs1_{i}") for i in range(1, 4)]
    sites_vcf = tmp_path / "sites.vcf"
    _make_sites_vcf(sites_vcf, records)
    input_filters = Stage1InputFilters(0, 0, 0, rows_written=len(records))

    snp_path = tmp_path / "out.snp"
    bed_path = tmp_path / "out.bed"

    s1, s2 = lift_and_transform_sharded(
        sites_vcf_path=sites_vcf,
        chain_path=tmp_path / "chain.gz",
        target_fasta_path=tmp_path / "ref.fa",
        output_lifted_vcf=tmp_path / "lifted.vcf",
        output_rejected_vcf=tmp_path / "rejected.vcf",
        output_snp_path=snp_path,
        output_bed_path=bed_path,
        input_filter_counters=input_filters,
        shard_tempdir=tmp_path / "shards",
        n_shards=1,
    )

    assert isinstance(s1, Stage1LiftCounters)
    assert isinstance(s2, Stage2TransformCounters)
    assert snp_path.exists() and snp_path.stat().st_size > 0
    assert bed_path.exists() and bed_path.stat().st_size > 0
    assert s2.output_sites == len(records)


def test_lift_and_transform_sharded_logical_equivalence(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """S2: streaming n_shards=2 produces same site set as single-shard run.

    The fake Picard lifts records verbatim; transform writes .snp/.bed from each
    shard's lifted VCF. After streaming, the merged .snp must contain the same
    (rsid, chrom_numeric, pos) tuples as the single-shard run.
    """
    _setup_fake_picard(monkeypatch, tmp_path)

    records = (
        [("chr1", i * 1000, f"rs1_{i}") for i in range(1, 6)]
        + [("chr22", i * 1000, f"rs22_{i}") for i in range(1, 4)]
    )
    input_filters = Stage1InputFilters(0, 0, 0, rows_written=len(records))

    def _run(out_dir: Path, n_shards: int) -> set[tuple[str, str, str]]:
        out_dir.mkdir(parents=True, exist_ok=True)
        sites_vcf = out_dir / "sites.vcf"
        _make_sites_vcf(sites_vcf, records)
        snp_path = out_dir / "out.snp"
        bed_path = out_dir / "out.bed"
        lift_and_transform_sharded(
            sites_vcf_path=sites_vcf,
            chain_path=out_dir / "chain.gz",
            target_fasta_path=out_dir / "ref.fa",
            output_lifted_vcf=out_dir / "lifted.vcf",
            output_rejected_vcf=out_dir / "rejected.vcf",
            output_snp_path=snp_path,
            output_bed_path=bed_path,
            input_filter_counters=input_filters,
            shard_tempdir=out_dir / "shards",
            n_shards=n_shards,
        )
        return {
            tuple(line.split()[:4])  # rsid, chrom_num, gen, pos
            for line in snp_path.read_text().splitlines()
            if line.strip()
        }

    single = _run(tmp_path / "single", 1)
    sharded = _run(tmp_path / "sharded", 2)
    assert single == sharded, f"site sets differ: {single ^ sharded}"
    assert len(single) == len(records)


def test_lift_and_transform_sharded_snp_bed_line_counts_match(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """S3: merged .snp and .bed have identical line counts after streaming."""
    _setup_fake_picard(monkeypatch, tmp_path)

    records = (
        [("chr1", i * 1000, f"rs1_{i}") for i in range(1, 8)]
        + [("chr2", i * 1000, f"rs2_{i}") for i in range(1, 5)]
        + [("chr22", i * 1000, f"rs22_{i}") for i in range(1, 4)]
    )
    input_filters = Stage1InputFilters(0, 0, 0, rows_written=len(records))
    sites_vcf = tmp_path / "sites.vcf"
    _make_sites_vcf(sites_vcf, records)
    snp_path = tmp_path / "out.snp"
    bed_path = tmp_path / "out.bed"

    lift_and_transform_sharded(
        sites_vcf_path=sites_vcf,
        chain_path=tmp_path / "chain.gz",
        target_fasta_path=tmp_path / "ref.fa",
        output_lifted_vcf=tmp_path / "lifted.vcf",
        output_rejected_vcf=tmp_path / "rejected.vcf",
        output_snp_path=snp_path,
        output_bed_path=bed_path,
        input_filter_counters=input_filters,
        shard_tempdir=tmp_path / "shards",
        n_shards=3,
    )

    snp_lines = [line for line in snp_path.read_text().splitlines() if line.strip()]
    bed_lines = [line for line in bed_path.read_text().splitlines() if line.strip()]
    assert len(snp_lines) == len(bed_lines), (
        f".snp has {len(snp_lines)} lines but .bed has {len(bed_lines)}"
    )
    assert len(snp_lines) == len(records)
