# Changelog

All notable changes to pileup-aadr will be documented in this file. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added (Day 1 — 2026-05-12)
- Initial repository scaffold per LLD §2 (repo skeleton)
- `pyproject.toml` with build-system + dependencies + dev extras + ruff/mypy/pytest config
- `pileup_aadr/counters.py` — 7 counter dataclasses (Stage1LiftCounters, Stage2TransformCounters, PileupCallerSummary, Stage3CallCounters, Stage4RejoinCounters, CoverageCounters, ExtractCounters)
- `pileup_aadr/errors.py` — `PileupAadrError` base + 20 named error classes per HLD §"Error class taxonomy"
- `pileup_aadr/logging_config.py` — stdlib logging setup + `JsonLinesFormatter`
- `pileup_aadr/format_detect.py` — BAM/CRAM detection, hg19/hg38 build detection, sample-name resolution, AADR `.snp` parser with duplicate-rsID invariant
- `pileup_aadr/chrom_lengths.py` — HG19/HG38 chromosome-length tables + CHROM_ORDER
- `pileup_aadr/cli.py` — click root group + `main()` entry point + error formatter
- `pileup_aadr/extract_cmd.py` — click decorator stub (orchestrator deferred to Days 3–5)
- `pileup_aadr/validate_cmd.py` + `validate_impl.py` — full 10-check pre-flight implementation
- `pileup_aadr/inspect_cmd.py` + `inspect_impl.py` — pure-Python AADR `.snp` summary
- `tests/conftest.py` + smoke tests for the implemented modules

### Added (Day 2 — 2026-05-12)
- Bundled UCSC `hg19ToHg38.over.chain.gz` (227,698 bytes) + SHA-256 sidecar
  at `pileup_aadr/data/`. Real chain bytes shipped; SHA pinned.
- `tool_wrapper.py` — `ToolSpec` + `ToolRunResult` + `ToolWrapper` class with
  binary lookup, version probe + caching, version comparison, JVM-args
  injection for jar specs (H9 fix), `run()` with stderr-to-disk discipline,
  `pipe()` for Stage 3 mpileup→pileupCaller chain with SIGPIPE handling.
  Five constants: `SAMTOOLS_SPEC`, `PILEUPCALLER_SPEC`, `PICARD_SPEC`,
  `JAVA_SPEC`, `MOSDEPTH_SPEC`. Picard JAR resolved via `_resolve_picard_jar()`
  with conda paths first (B1 fix).
- `lift.py` (partial — chain resolution only; Stage 1 lift implementation
  lands Day 3) — `get_bundled_chain_path()` with always-on SHA verification,
  3-tier `chain_file_path()` resolver (--chain → env → bundled),
  `resolve_chain_for_extract()` env-aware wrapper, `_verify_user_chain_sha()`
  for `--strict-chain-sha` enforcement.
- `ref_resolve.py` — `resolve_ref_fasta()` 3-tier resolver
  (--ref-fasta → env → BAM @PG), `_extract_ref_from_bam_pg()` regex-based
  extraction from BAM header lines, `verify_fasta_matches_bam_build()` with
  closest-match logic (hg19 and hg38 chr1 differ by only 294 KB so simple
  ±1 Mb tolerance overlaps; closest-match disambiguates correctly).
- `validate_impl.py` updates — Day-1 binary-presence stubs replaced with
  real version probes via `ToolWrapper`. Chain check uses `resolve_chain_for_extract`
  (verifies bundled SHA). Ref FASTA check uses `verify_fasta_matches_bam_build`
  for proper build-mismatch diagnosis.

### Tests added (Day 2; 120 total now)
- `test_lift_chain.py` — 11 tests covering bundled SHA matches sidecar,
  canonical-size sanity check (200-250 KB), 3-tier resolution matrix,
  --strict-chain-sha + --insecure-chain interactions, env-var resolver.
- `test_ref_resolve.py` — 13 tests covering `verify_fasta_matches_bam_build`
  (hg19/hg38 OK, build mismatch, missing .fai, no chr1 in .fai, within-tolerance
  drift), `_extract_ref_from_bam_pg` (success, no @PG, stale path), full
  3-tier `resolve_ref_fasta` flow.
- `test_tool_wrapper.py` — 16 tests covering all 5 spec constants,
  binary-not-found errors per tool, --version regex parsing, version-too-old
  rejection, version-cache behavior, `_build_invocation` for jar + non-jar
  + jvm_args (H9), ToolRunResult/ToolSpec frozen-ness.

### Bugs surfaced + fixed during testing
- `verify_fasta_matches_bam_build` used the same overlap-tolerance bug as
  Day 1's `detect_aadr_build` (hg19/hg38 chr1 are 294 KB apart, well within
  ±1 Mb). Switched to closest-match logic to disambiguate.
- Test fixture for "patch-level drift" used `HG38_CHR1_LENGTH + 500_000`
  which actually crosses into hg19 territory (500K > 294K gap). Updated to
  `HG38_CHR1_LENGTH - 100_000` which stays unambiguously closer to hg38.
- `ToolWrapper._version_cache` triggered `RUF012` (mutable class attr);
  annotated with `ClassVar` to make the intent explicit.
