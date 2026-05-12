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

### Outstanding
- Stage 1–4 implementations (Days 2–5)
- Synthetic test fixtures (BAM, AADR slice, synthetic chains; per LLD §18)
- Bundled chain file SHA pinning (placeholder for Day 1)
