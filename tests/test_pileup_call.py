"""Tests for pileup_call.py — Stage 3 (mpileup | pileupCaller pipe).

The Picard-equivalent integration tests mock `ToolWrapper.pipe` to return canned
exit codes + write a captured pileupCaller stderr fixture. Tests cover: clean
run → counters populated; downstream non-zero → ToolSubprocessError; upstream
non-zero (other than SIGPIPE-141) → ToolSubprocessError; upstream 141 +
downstream 0 → tolerated; stderr-parser format paths.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from pileup_aadr.counters import PileupCallerSummary, Stage3CallCounters
from pileup_aadr.errors import PileupAadrInternalError, ToolSubprocessError
from pileup_aadr.pileup_call import (
    parse_pileupcaller_stderr,
    run_pileup_call,
    run_pileup_call_shards,
)
from pileup_aadr.tool_wrapper import ToolRunResult

FIXTURES_STDERR = Path(__file__).parent / "fixtures" / "stderr"


# --- parse_pileupcaller_stderr ---


def test_parse_clean_run() -> None:
    """Captured 1131600-site pileupCaller stderr parses to PileupCallerSummary."""
    text = (FIXTURES_STDERR / "pileupcaller_clean.stderr").read_text()
    summary = parse_pileupcaller_stderr(text)
    assert summary.total_sites == 1131600
    assert summary.non_missing_calls == 1107213
    assert summary.avg_raw_reads == 21.4
    assert summary.avg_damage_cleaned_reads == 21.4
    assert summary.avg_sampled_from == 21.4


def test_parse_scientific_notation_in_avg_columns() -> None:
    """Captured low-coverage stderr with avgRawReads in scientific notation
    (e.g., 3.68e-2) must parse — bug-fix regression for the smoke test on
    the tiny.bam where coverage was tiny enough to trigger e-notation."""
    text = (FIXTURES_STDERR / "pileupcaller_scientific.stderr").read_text()
    summary = parse_pileupcaller_stderr(text)
    assert summary.total_sites == 4965
    assert summary.non_missing_calls == 9
    assert summary.avg_raw_reads == 0.036858006042296075


def test_parse_missing_header_raises() -> None:
    """Stderr lacking the SampleName header → PileupAadrInternalError."""
    text = "Some unrelated output without the header.\n"
    with pytest.raises(PileupAadrInternalError, match="header"):
        parse_pileupcaller_stderr(text)


def test_parse_missing_data_line_raises() -> None:
    """Header present but no data line → PileupAadrInternalError."""
    text = (
        "SampleName\tTotalSites\tNonMissingCalls\tavgRawReads"
        "\tavgDamageCleanedReads\tavgSampledFrom\n"
        "# (no data row at all)\n"
    )
    with pytest.raises(PileupAadrInternalError, match="data line"):
        parse_pileupcaller_stderr(text)


# --- run_pileup_call (pipe mocked) ---


def _setup_run_mocks(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    *,
    upstream_exit: int,
    downstream_exit: int,
    pileupcaller_stderr_text: str,
    mpileup_stderr_text: str = "",
) -> None:
    """Patch subprocess.run (version probe) + ToolWrapper.pipe + _check_version."""
    import subprocess
    from unittest.mock import MagicMock

    # samtools/pileupCaller version probes both stub out
    monkeypatch.setattr(
        subprocess, "run",
        lambda *_a, **_kw: MagicMock(stdout="1.6.0.0", stderr="", returncode=0),
    )

    from pileup_aadr import pileup_call

    # Bypass binary-on-PATH lookup — neither samtools nor pileupCaller may be in
    # the test env. _resolve_binary returns a Path; the actual binary is never
    # invoked because we replace `pipe()` below.
    monkeypatch.setattr(
        pileup_call.ToolWrapper, "_resolve_binary",
        lambda _self, spec: Path(f"/usr/bin/fake_{spec.binary}"),
    )

    def fake_pipe(
        self: object,
        downstream: object,
        *,
        upstream_args: list[str],
        downstream_args: list[str],
        upstream_stderr_to: Path,
        downstream_stderr_to: Path,
    ) -> tuple[ToolRunResult, ToolRunResult]:
        upstream_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        upstream_stderr_to.write_text(mpileup_stderr_text)
        downstream_stderr_to.write_text(pileupcaller_stderr_text)
        return (
            ToolRunResult(
                exit_code=upstream_exit,
                stdout=None,
                stderr_path=upstream_stderr_to,
                stderr_text=None,
                wallclock_seconds=0.1,
                peak_rss_mb=None,
            ),
            ToolRunResult(
                exit_code=downstream_exit,
                stdout=None,
                stderr_path=downstream_stderr_to,
                stderr_text=None,
                wallclock_seconds=0.1,
                peak_rss_mb=None,
            ),
        )

    monkeypatch.setattr(pileup_call.ToolWrapper, "pipe", fake_pipe)
    monkeypatch.setattr(pileup_call.ToolWrapper, "_check_version", lambda _self: None)


def _common_run_args(tmp_path: Path) -> dict[str, Any]:
    bam = tmp_path / "user.bam"
    bam.touch()
    snp = tmp_path / "aadr.snp"
    snp.touch()
    bed = tmp_path / "aadr.bed"
    bed.touch()
    fasta = tmp_path / "ref.fa"
    fasta.touch()
    return {
        "bam_path": bam,
        "snp_path": snp,
        "bed_path": bed,
        "target_fasta_path": fasta,
        "output_prefix": tmp_path / "call" / "user_hg38",
        "sample_name": "Carsten",
        "pop_name": "TestPop",
    }


def test_clean_run_returns_counters(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Both processes exit 0 + pileupCaller stderr parses → Stage3CallCounters."""
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=0, downstream_exit=0,
        pileupcaller_stderr_text=(
            FIXTURES_STDERR / "pileupcaller_clean.stderr"
        ).read_text(),
    )
    counters = run_pileup_call(**_common_run_args(tmp_path))
    assert counters.pileupcaller_summary.total_sites == 1131600
    assert counters.pileupcaller_summary.non_missing_calls == 1107213
    assert counters.wallclock_seconds >= 0


def test_downstream_nonzero_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """pileupCaller exit != 0 → ToolSubprocessError naming pileupCaller."""
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=0, downstream_exit=1,
        pileupcaller_stderr_text="ERROR: bad input file\n",
    )
    with pytest.raises(ToolSubprocessError) as excinfo:
        run_pileup_call(**_common_run_args(tmp_path))
    assert "pileupCaller" in excinfo.value.what
    assert "exit code 1" in excinfo.value.why


def test_upstream_sigpipe_with_clean_downstream_tolerated(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Upstream exits 141 (SIGPIPE) + downstream exits 0 → tolerated, no raise."""
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=141, downstream_exit=0,
        pileupcaller_stderr_text=(
            FIXTURES_STDERR / "pileupcaller_clean.stderr"
        ).read_text(),
    )
    counters = run_pileup_call(**_common_run_args(tmp_path))
    assert counters.pileupcaller_summary.total_sites == 1131600


def test_upstream_real_error_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Upstream exit != 0 and != 141 → ToolSubprocessError naming samtools."""
    _setup_run_mocks(
        monkeypatch, tmp_path,
        upstream_exit=1, downstream_exit=0,
        pileupcaller_stderr_text=(
            FIXTURES_STDERR / "pileupcaller_clean.stderr"
        ).read_text(),
        mpileup_stderr_text="samtools: missing BAM index\n",
    )
    with pytest.raises(ToolSubprocessError) as excinfo:
        run_pileup_call(**_common_run_args(tmp_path))
    assert "samtools mpileup" in excinfo.value.what
    assert "missing BAM index" in excinfo.value.why


def test_samtools_args_never_include_dash_at(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Issue #2 regression: the constructed `samtools mpileup` argv MUST NOT
    contain `-@` (mpileup rejects it; v0.1.0-0.1.1 mistakenly passed it)."""
    captured: dict[str, list[str]] = {}

    def capture_pipe(
        self: object,
        downstream: object,
        *,
        upstream_args: list[str],
        downstream_args: list[str],
        upstream_stderr_to: Path,
        downstream_stderr_to: Path,
    ) -> tuple[ToolRunResult, ToolRunResult]:
        captured["upstream_args"] = upstream_args
        upstream_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        upstream_stderr_to.write_text("")
        downstream_stderr_to.write_text(
            (FIXTURES_STDERR / "pileupcaller_clean.stderr").read_text(),
        )
        return (
            ToolRunResult(exit_code=0, stdout=None, stderr_path=upstream_stderr_to,
                          stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None),
            ToolRunResult(exit_code=0, stdout=None, stderr_path=downstream_stderr_to,
                          stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None),
        )

    from pileup_aadr import pileup_call
    monkeypatch.setattr(
        pileup_call.ToolWrapper, "_resolve_binary",
        lambda _self, spec: Path(f"/usr/bin/fake_{spec.binary}"),
    )
    monkeypatch.setattr(pileup_call.ToolWrapper, "_check_version", lambda _self: None)
    monkeypatch.setattr(pileup_call.ToolWrapper, "pipe", capture_pipe)

    run_pileup_call(**_common_run_args(tmp_path))
    assert "-@" not in captured["upstream_args"]


def test_region_argument_appended_to_samtools_args(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """region='chr1' appears in the samtools mpileup argv (shard-scoped call)."""
    captured: dict[str, list[str]] = {}

    def capture_pipe(
        self: object,
        downstream: object,
        *,
        upstream_args: list[str],
        downstream_args: list[str],
        upstream_stderr_to: Path,
        downstream_stderr_to: Path,
    ) -> tuple[ToolRunResult, ToolRunResult]:
        captured["upstream_args"] = upstream_args
        upstream_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        upstream_stderr_to.write_text("")
        downstream_stderr_to.write_text(
            (FIXTURES_STDERR / "pileupcaller_clean.stderr").read_text()
        )
        return (
            ToolRunResult(exit_code=0, stdout=None, stderr_path=upstream_stderr_to,
                          stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None),
            ToolRunResult(exit_code=0, stdout=None, stderr_path=downstream_stderr_to,
                          stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None),
        )

    from pileup_aadr import pileup_call as pc_mod
    monkeypatch.setattr(pc_mod.ToolWrapper, "_resolve_binary",
                        lambda _self, spec: Path(f"/usr/bin/fake_{spec.binary}"))
    monkeypatch.setattr(pc_mod.ToolWrapper, "_check_version", lambda _self: None)
    monkeypatch.setattr(pc_mod.ToolWrapper, "pipe", capture_pipe)

    run_pileup_call(**_common_run_args(tmp_path), region="chr1")
    assert "chr1" in captured["upstream_args"]


def test_region_comes_after_bam_in_samtools_args(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """samtools mpileup requires the region AFTER the BAM file: `mpileup [opts] bam [region]`.
    Region before BAM causes samtools to try to open the region string as a BAM file."""
    captured: dict[str, list[str]] = {}

    def capture_pipe(
        self: object,
        downstream: object,
        *,
        upstream_args: list[str],
        downstream_args: list[str],
        upstream_stderr_to: Path,
        downstream_stderr_to: Path,
    ) -> tuple[ToolRunResult, ToolRunResult]:
        captured["upstream_args"] = upstream_args
        upstream_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        upstream_stderr_to.write_text("")
        downstream_stderr_to.write_text(
            (FIXTURES_STDERR / "pileupcaller_clean.stderr").read_text()
        )
        return (
            ToolRunResult(exit_code=0, stdout=None, stderr_path=upstream_stderr_to,
                          stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None),
            ToolRunResult(exit_code=0, stdout=None, stderr_path=downstream_stderr_to,
                          stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None),
        )

    from pileup_aadr import pileup_call as pc_mod
    monkeypatch.setattr(pc_mod.ToolWrapper, "_resolve_binary",
                        lambda _self, spec: Path(f"/usr/bin/fake_{spec.binary}"))
    monkeypatch.setattr(pc_mod.ToolWrapper, "_check_version", lambda _self: None)
    monkeypatch.setattr(pc_mod.ToolWrapper, "pipe", capture_pipe)

    args = _common_run_args(tmp_path)
    run_pileup_call(**args, region="chr1")

    argv = captured["upstream_args"]
    bam_str = str(args["bam_path"])
    assert bam_str in argv
    assert "chr1" in argv
    assert argv.index(bam_str) < argv.index("chr1"), (
        "region must come AFTER BAM in samtools mpileup argv; "
        f"got: {argv}"
    )


# --- calling mode → pileupCaller flag ---


def _capture_downstream_args(monkeypatch: pytest.MonkeyPatch) -> dict[str, list[str]]:
    """Patch ToolWrapper so the pileupCaller (downstream) argv is captured."""
    captured: dict[str, list[str]] = {}

    def capture_pipe(
        self: object,
        downstream: object,
        *,
        upstream_args: list[str],
        downstream_args: list[str],
        upstream_stderr_to: Path,
        downstream_stderr_to: Path,
    ) -> tuple[ToolRunResult, ToolRunResult]:
        captured["downstream_args"] = downstream_args
        upstream_stderr_to.parent.mkdir(parents=True, exist_ok=True)
        upstream_stderr_to.write_text("")
        downstream_stderr_to.write_text(
            (FIXTURES_STDERR / "pileupcaller_clean.stderr").read_text()
        )
        return (
            ToolRunResult(exit_code=0, stdout=None, stderr_path=upstream_stderr_to,
                          stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None),
            ToolRunResult(exit_code=0, stdout=None, stderr_path=downstream_stderr_to,
                          stderr_text=None, wallclock_seconds=0.1, peak_rss_mb=None),
        )

    from pileup_aadr import pileup_call as pc_mod
    monkeypatch.setattr(pc_mod.ToolWrapper, "_resolve_binary",
                        lambda _self, spec: Path(f"/usr/bin/fake_{spec.binary}"))
    monkeypatch.setattr(pc_mod.ToolWrapper, "_check_version", lambda _self: None)
    monkeypatch.setattr(pc_mod.ToolWrapper, "pipe", capture_pipe)
    return captured


def test_default_calling_mode_is_random_haploid(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """No calling_mode → pileupCaller gets --randomHaploid (matches the AADR panel)."""
    captured = _capture_downstream_args(monkeypatch)
    run_pileup_call(**_common_run_args(tmp_path))
    argv = captured["downstream_args"]
    assert "--randomHaploid" in argv
    assert "--randomDiploid" not in argv
    assert "--majorityCall" not in argv
    assert "--seed" in argv


@pytest.mark.parametrize("mode", ["randomHaploid", "randomDiploid", "majorityCall"])
def test_calling_mode_maps_to_pileupcaller_flag(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, mode: str,
) -> None:
    """calling_mode selects the matching --<mode> flag; --seed is passed for every
    mode (majorityCall uses it to break equal-depth ties → reproducibility)."""
    captured = _capture_downstream_args(monkeypatch)
    run_pileup_call(**_common_run_args(tmp_path), calling_mode=mode)
    argv = captured["downstream_args"]
    assert f"--{mode}" in argv
    other_modes = {"randomHaploid", "randomDiploid", "majorityCall"} - {mode}
    for other in other_modes:
        assert f"--{other}" not in argv
    assert "--seed" in argv


# --- run_pileup_call_shards ---


def _write_two_chrom_sites(snp_path: Path, bed_path: Path) -> None:
    snp_path.parent.mkdir(parents=True, exist_ok=True)
    snp_path.write_text("rs1\t1\t0.0\t1000\tA\tG\nrs22\t22\t0.0\t2000\tC\tT\n")
    bed_path.write_text("chr1\t999\t1000\nchr22\t1999\t2000\n")


def _make_fake_run_pileup_call(monkeypatch: pytest.MonkeyPatch) -> list[dict]:
    """Patch module-level run_pileup_call to write minimal shard outputs; return call log."""
    from pileup_aadr import pileup_call as pc_mod

    calls: list[dict] = []

    def fake_run(**kwargs):  # type: ignore[no-untyped-def]
        calls.append(dict(kwargs))
        snp_lines = kwargs["snp_path"].read_text().splitlines()
        op = kwargs["output_prefix"]
        op.parent.mkdir(parents=True, exist_ok=True)
        with open(f"{op}.geno", "w") as gh, open(f"{op}.snp", "w") as sh:
            for line in snp_lines:
                gh.write("0\n")
                sh.write(line + "\n")
        Path(f"{op}.ind").write_text(f"{kwargs['sample_name']}\tU\t{kwargs['pop_name']}\n")
        return Stage3CallCounters(
            wallclock_seconds=0.1,
            pileupcaller_summary=PileupCallerSummary(
                total_sites=len(snp_lines), non_missing_calls=len(snp_lines),
                avg_raw_reads=5.0, avg_damage_cleaned_reads=5.0, avg_sampled_from=5.0,
            ),
        )

    monkeypatch.setattr(pc_mod, "run_pileup_call", fake_run)
    return calls


def test_shards_threads1_single_process_no_fanout(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """threads=1 → single run_pileup_call invocation, per_shard=[]."""
    calls = _make_fake_run_pileup_call(monkeypatch)

    snp = tmp_path / "sites.snp"
    bed = tmp_path / "sites.bed"
    _write_two_chrom_sites(snp, bed)

    counters = run_pileup_call_shards(
        bam_path=tmp_path / "user.bam",
        sites_snp_path=snp, sites_bed_path=bed,
        target_fasta_path=tmp_path / "ref.fa",
        output_prefix=tmp_path / "call" / "out",
        sample_name="S", pop_name="P",
        shard_dir=tmp_path / "shards",
        master_seed=42, threads=1,
    )

    assert len(calls) == 1
    assert counters.per_shard == []
    assert calls[0]["seed"] == 42
    assert calls[0].get("region") is None


def test_shards_threads2_fans_out_and_merges(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """threads=2 with 2 chromosomes → 2 per_shard entries, merged output."""
    _make_fake_run_pileup_call(monkeypatch)

    snp = tmp_path / "sites.snp"
    bed = tmp_path / "sites.bed"
    _write_two_chrom_sites(snp, bed)

    counters = run_pileup_call_shards(
        bam_path=tmp_path / "user.bam",
        sites_snp_path=snp, sites_bed_path=bed,
        target_fasta_path=tmp_path / "ref.fa",
        output_prefix=tmp_path / "out",
        sample_name="S", pop_name="P",
        shard_dir=tmp_path / "shards",
        master_seed=42, threads=2,
    )

    assert len(counters.per_shard) == 2
    shard_chroms = {s.chromosome for s in counters.per_shard}
    assert shard_chroms == {"chr1", "chr22"}
    geno_lines = Path(f"{tmp_path / 'out'}.geno").read_text().splitlines()
    assert len(geno_lines) == 2  # one per site


def test_shards_failure_propagates(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Shard run_pileup_call failure → ToolSubprocessError propagated."""
    from pileup_aadr import pileup_call as pc_mod

    def fail_run(**_kw):  # type: ignore[no-untyped-def]
        raise ToolSubprocessError(what="pileupCaller", why="shard crashed", fix="check logs")

    monkeypatch.setattr(pc_mod, "run_pileup_call", fail_run)

    snp = tmp_path / "sites.snp"
    bed = tmp_path / "sites.bed"
    snp.parent.mkdir(parents=True, exist_ok=True)
    snp.write_text("rs1\t1\t0.0\t1000\tA\tG\n")
    bed.write_text("chr1\t999\t1000\n")

    with pytest.raises(ToolSubprocessError, match="shard crashed"):
        run_pileup_call_shards(
            bam_path=tmp_path / "user.bam",
            sites_snp_path=snp, sites_bed_path=bed,
            target_fasta_path=tmp_path / "ref.fa",
            output_prefix=tmp_path / "out",
            sample_name="S", pop_name="P",
            shard_dir=tmp_path / "shards",
            master_seed=42, threads=2,
        )


def test_shards_per_shard_region_matches_chromosome(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Each per-chromosome shard is called with its chromosome as the region arg."""
    calls = _make_fake_run_pileup_call(monkeypatch)

    snp = tmp_path / "sites.snp"
    bed = tmp_path / "sites.bed"
    _write_two_chrom_sites(snp, bed)

    run_pileup_call_shards(
        bam_path=tmp_path / "user.bam",
        sites_snp_path=snp, sites_bed_path=bed,
        target_fasta_path=tmp_path / "ref.fa",
        output_prefix=tmp_path / "out",
        sample_name="S", pop_name="P",
        shard_dir=tmp_path / "shards",
        master_seed=42, threads=2,
    )

    regions = {c.get("region") for c in calls}
    assert regions == {"chr1", "chr22"}
