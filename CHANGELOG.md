# Changelog

All notable changes to pileup-aadr will be documented in this file. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] — 2026-05-12

Beta promotion. Six of the seven v0.2 items from the post-v0.1
reconciliation review landed; the seventh (region-parallel mpileup
fan-out) is the largest single item and rolls forward to v0.3 along
with any other Stage-3 wallclock improvements that profiling surfaces.
The headline change vs v0.1: pileup-aadr is now genuinely panel-
agnostic — HumanOrigins, chrY-only haplogroup, chrM-only mtDNA, and
arbitrary custom EIGENSTRAT panels all work end-to-end without
crashing on the 1240k-specific autosomal coverage gate. The HLD
spec'd LLD-#19 full-chain f2 equivalence test is also wired (closes
"verified equivalent to mergeit's pipeline").

249 tests passing across the 6-cell + bio-tools CI matrix; ruff
clean. See per-day Added/Fixed sections below for the full v0.2
implementation history.

### Added (v0.2 in progress)
- **Panel classification + non-autosomal coverage-gate skip** (HLD §"Out-of-
  scope reachable extensions": HO panel + chrY-only + chrM-only + custom
  panels). New `format_detect.classify_aadr_chrom_set(df) -> str` tags the
  panel as one of `autosomes+sex` / `autosomes_only` / `sex_only` /
  `chrY_only` / `chrM_only` / `custom`. The orchestrator uses this in
  `_evaluate_coverage_gate` to skip the autosomal threshold cleanly with
  an INFO log when the panel is non-autosomal — chrY haplogroup and
  chrM mtDNA workflows no longer fail with a misleading "below
  --min-coverage 500000" message. Per-chrom counts still flow into
  `coverage.per_chrom_call_count` for downstream sanity. `inspect` adds
  `chrom_set` to its summary fields. README documents the panel classes
  and their gate semantics. 7 new tests (6 classifier + 1 orchestrator
  gate-skip end-to-end against a chrY-only synthetic panel).

  Practical effect: the tool was already panel-agnostic by construction
  (any EIGENSTRAT `.snp` reads cleanly), but the autosomal coverage gate
  was a 1240k-specific assumption baked into the orchestrator. Removing
  that assumption opens HumanOrigins (~600K sites; user picks own
  threshold), chrY-only haplogroup workflows, chrM-only mtDNA workflows,
  and arbitrary custom panels (with documentation).

- **chr20-as-anchor fallback for build detection** (customer issue #1
  follow-up suggestion). When chr1 is absent from the BAM `@SQ` headers
  or AADR `.snp` (e.g., chrY-only haplogroup workflows or a chr-by-chr
  AADR slice), `detect_bam_build` and `detect_aadr_build` now fall back
  to chr20 (hg19=63,025,520 / hg38=64,444,167; the wider 1.4 Mb gap
  disambiguates cleanly via the same closest-match logic). Diagnostic
  on total failure names BOTH anchors so the user sees what was tried.
  6 new tests in `test_format_detect.py` covering both anchors × both
  builds + the no-anchor failure path.
- **`validate` flag-probe checks** (customer issue #2 follow-up
  suggestion). New checks run each tool's `--help` and grep the output
  for the specific flags pileup-aadr will use:
    - `samtools mpileup`: `-B`, `-q`, `-Q`, `-R`, `-f`, `-l`
    - `pileupCaller`: `--randomDiploid`, `--seed`, `--sampleNames`,
      `--samplePopName`, `-e`
    - `mosdepth`: `--threads`, `--no-per-base`, `--quantize`, `--by`
    - `picard.jar LiftoverVcf` (skipped on no-lift fast path):
      `--INPUT`, `--OUTPUT`, `--CHAIN`, `--REFERENCE_SEQUENCE`,
      `--REJECT`, `--RECOVER_SWAPPED_REF_ALT`, `--MAX_RECORDS_IN_RAM`
  Catches the bug class where samtools/Picard rename or drop a flag in
  a new release. Word-boundary matching (so `-q` doesn't false-pass on
  `-quiet`). Missing binary or hung subprocess fails cleanly with a
  clear message rather than a traceback. 4 new unit tests in
  `test_validate_cmd.py` covering the pass path, missing-flag failure,
  word-boundary discipline, missing-binary failure.
- **LLD #19 full-chain f2 equivalence test.** Two layers gated on AT2 fork
  (`carstenerickson/admixtools@production/v1.0`) + pgen-samplebind:
    - **Layer A** (always runs when toolchain present, ~5s): loads a
      pre-extracted pileup-aadr triplet from `tests/integration/fixtures/
      lld19/pileup_aadr_output/`, merges with the AADR chr22 5-pop
      reference subset via `pgen-samplebind merge`, runs AT2 `extract_f2`
      on the merged PFILE, asserts the f2 array matches the frozen
      `f2_baseline.rds` within `max_dev < 1e-9`.
    - **Layer B** (gated on `PILEUP_AADR_LLD19_BAM` + `_REF` env vars,
      ~30s): full chain BAM → pileup-aadr extract → samplebind merge →
      AT2 extract_f2 → diff baseline. Catches drift in pileup-aadr
      itself.
- `tests/integration/fixtures/lld19/` — 1.8 MB committed fixtures: AADR
  chr22 5-pop subset (Loschbour, Stuttgart, MA1, Denisova, Altai
  Neanderthal; ~16K SNPs × 11 samples), a frozen pileup-aadr extract
  output for layer A, the 6×6×13 f2 baseline `.rds`, and the original
  extract `--report-json` for provenance. README documents the panel
  composition + sources.
- `scripts/regen_lld19_baseline.sh` — end-to-end script regenerating
  every fixture from a clean toolchain (AADR PUB → 5-pop subset → BAM
  subset → pileup-aadr extract → pgen-samplebind merge → AT2
  extract_f2). For when pileup-aadr's pipeline legitimately moves the
  numbers (new Stage-4 invariant, etc.).
- `tests/test_extract_cmd.py` — 3 click-layer contract tests:
  default-no-flag-disables-mpileup-BAQ, --enable-baq-enables-it,
  CLI-default-equals-dataclass-default. Catches the BAQ-flip class of
  bug surfaced below.

### Fixed (v0.2 in progress)
- **CLI `--enable-baq` flip was inverted.** v0.1.0-0.1.2's
  `extract_cmd.extract` did `kwargs["no_baq"] = not enable_baq`, which
  produced **BAQ-ENABLED** mpileup as the default — the opposite of
  HLD §"CLI reference > pileup / call" which specifies "default: -B is
  passed, disabling samtools BAQ to match pileupCaller's recommended
  cmdline". The dataclass default `no_baq=False` was correct, so direct
  programmatic instantiation produced the right behavior — but
  `pileup-aadr extract` from the CLI without `--enable-baq` produced a
  different `.geno` than `run_extract(ExtractCliArgs(...))` would. This
  divergence is what broke LLD #19 layer-B at first run. Fix:
  `kwargs["no_baq"] = enable_baq` (pass-through, not negation).
  Regression tests in `test_extract_cmd.py`.
- README documents the local-test toolchain setup
  (`PICARD_JAR=`, `openjdk` on PATH, `pgen-samplebind`'s venv on PATH,
  AT2 fork install) for running `tests/integration/` locally rather
  than only on CI's bio-tools job.

## [0.1.2] — 2026-05-12

Two real bug fixes from the pipeline customer.

### Fixed
- **#1 — `detect_bam_build` misclassified hg38 BAMs as hg19.** hg19 chr1
  (249,250,621) and hg38 chr1 (248,956,422) are only 294 KB apart — well
  INSIDE the ±1 Mb tolerance window. The first-match-wins logic always
  returned `"hg19"` for hg38 BAMs because the hg19 check fired first
  within tolerance. Switched to closest-match (same fix already applied
  to `detect_aadr_build` + `verify_fasta_matches_bam_build`; this was a
  regression of the same bug class). Workaround pre-fix: pass
  `--bam-build hg38` explicitly.
- **#2 — `samtools mpileup -@` is invalid.** `mpileup` is single-threaded
  and has no `-@`/`--threads` flag in any release through 1.23.1
  (verified: `samtools view`'s `--threads` is a separate codepath).
  v0.1.0-0.1.1 erroneously appended `-@ N` for `--threads > 1`
  invocations, causing samtools to error ~22s into Stage 3. The
  `--threads` CLI flag is preserved for back-compat but now treated as a
  no-op; pileupCaller is also single-threaded (Haskell) so Stage 3 is
  always single-threaded by tool nature. A WARN log fires when
  `--threads > 1` to make the no-op visible. The `--no-thread-cap` flag
  is similarly preserved as a no-op. Workaround pre-fix: pass
  `--threads 1`.

### Tests added
- `test_format_detect.test_detect_bam_build_hg38_not_misclassified_as_hg19`
  is the #1 regression test (asserts the closest-match logic returns
  `"hg38"` for chr1=248_956_422), plus 3 sibling tests covering hg19
  detection, override short-circuit, and unknown-assembly raises.
- `test_pileup_call.test_samtools_args_never_include_dash_at` is the #2
  regression test (asserts the constructed mpileup argv never contains
  `-@` for any `--threads` value: 1, 4, 16). Plus
  `test_threads_gt_1_warns_no_op` and `test_threads_eq_1_silent` cover
  the new WARN behavior.

## [0.1.1] — 2026-05-12

HLD/LLD reconciliation pass: closes the four functional gaps surfaced in
the v0.1.0 reconciliation review. Adds the bio-tools CI job that would
have caught the v0.1.0 smoke-test class of bug automatically. 224 tests
now (15 new), 6-cell + bio-tools CI matrix.

### Added
- `coverage` subcommand emits the spec'd 9 columns: `chrom`, `length`,
  `bases`, `mean_coverage`, `median_coverage`, `fraction_at_>=1x`,
  `fraction_at_>=5x`, `fraction_at_>=10x`, `fraction_at_>=30x`. Median
  + the four fraction columns are derived from
  `<prefix>.mosdepth.global.dist.txt` (cumulative depth distribution).
  Defensive: missing dist file → quantile cols set to NaN, no crash.
  `--quantize` is forwarded to mosdepth so `.quantized.bed.gz` is also
  written for users who want bin-level analysis. Closes HLD §"CLI
  reference > coverage" output-columns spec.
- `pileup_aadr/dict_resolve.py` — `.dict` auto-generation per HLD §"Chain
  & reference dependencies > Target FASTA sequence dictionary":
    - `find_existing_dict` checks both `<fasta>.dict` (Picard's default
      output convention) and `<fasta-stem>.dict` (GATK resource bundle).
    - `find_or_user_cache_dict_path` decides where to write a NEW `.dict`:
      alongside the FASTA when its parent dir is writable, else
      `~/.cache/pileup-aadr/dicts/<sha-of-abs-path>.dict`.
    - `ensure_target_fasta_dict` is the front door — returns the existing
      `.dict` if any, else generates via `picard CreateSequenceDictionary`
      + caches (~23s on hg38).
  Wired into `extract_orch.run_extract` after `resolve_ref_fasta` (skipped
  on the no-lift fast path; Picard isn't invoked there) + into
  `validate_impl` as a new `target FASTA .dict` check (PASS if dict exists,
  WARN if extract would auto-generate, SKIP on no-lift).
- `tests/integration/` — real-binary integration suite (skipped silently
  when binaries are missing). Catches the smoke-test class of bug
  (Picard --version regex, CompletedProcess.pid, scientific-notation
  parser drift) automatically rather than requiring manual GCP-instance
  re-discovery.
    - `test_real_binaries.py`: 4 version-probe tests + 1 end-to-end
      `extract` against UCSC hg38 chr22 (downloaded + cached) + 1
      `coverage`-via-real-mosdepth test.
    - `test_pgen_samplebind_composition.py`: 3 always-runs shape-contract
      tests verifying `<prefix>.pseudohaploid.json` carries every field
      pgen-samplebind#2 enumerates + 1 live-binary test (skipped without
      pgen-samplebind on PATH). Maps to LLD #19.
- `benchmarks/bench_extract.py` — reproducible end-to-end perf bench per
  LLD §"Performance benchmark (test #20)". Points at user-supplied bench
  data via `PILEUP_AADR_BENCH_BAM` / `_SNP` / `_REF` env vars; skipped
  when unset. Reports wallclock; asserts a generous 10-min ceiling
  (tighter regression detection is hardware-dependent, deferred to user
  pinning).
- `.github/workflows/ci.yml` — new `bio-tools` job: installs Picard 3.3.0
  + samtools 1.23.1 + sequencetools 1.6.0 + mosdepth 0.3.6 from bioconda
  via `conda-incubator/setup-miniconda@v3`, caches the chr22 FASTA
  fixture, runs `pytest tests/integration/`. Mirrors LLD §18's `current`
  matrix entry.
- `pytest` config: `python_files = ["test_*.py", "bench_*.py"]` so
  `pytest benchmarks/` discovers the bench module without treating bench
  files as part of the default unit suite.

### Fixed
- HLD CLI ref drift: `--no-baq` (HLD spec) → `--enable-baq` (impl). The
  positive-flag form is what shipped at v0.1.0; HLD now matches. The
  orchestrator + dataclass field stay positively-named (`no_baq`); the
  click layer flips at the boundary.

## [0.1.0] — 2026-05-12

First tagged release. Surface complete: 4 subcommands (`extract`,
`validate`, `coverage`, `inspect`) functional end-to-end; 209 unit +
integration tests passing on a 6-cell CI matrix (Python 3.11/3.12/3.13
× ubuntu-latest/macos-latest); ruff clean. Real-binary smoke test on
GCP n2-standard-16 verified the v2.1 captured baselines (Picard 3.1.1
+ samtools 1.23.1 + pileupCaller 1.6.0.0 against a 5000-site AADR
slice + hg38 BAM): 4965/4967 lifted (99.96% yield, 16 SwappedAlleles,
2 NoTarget) in 13.6s wallclock; the load-bearing
`Stage1.swapped_alleles_count == Stage4.ref_alt_swap_count` invariant
holds. See per-day Added/Fixed sections below for the per-stage
implementation history.

### Fixed (Day 8 smoke test — 2026-05-12)
Three real bugs surfaced by running the full pipeline end-to-end against
real Picard 3.1.1 / samtools 1.23.1 / pileupCaller 1.6.0.0 / hg38 FASTA on
a 5000-site AADR slice + an hg38 BAM. The mocked-subprocess test suite
didn't catch any of them; the smoke test was the only path that could.

- `tool_wrapper.PICARD_SPEC.version_args` was `["--version"]`, but
  `java -jar picard.jar --version` (no subcommand) prints help + exits 1
  — the version string `Version:X.Y.Z` only appears when `--version`
  comes AFTER a subcommand. Switched to `["LiftoverVcf", "--version"]`
  since LiftoverVcf is what we actually invoke. Validate's picard.jar
  check now PASSes against Picard 3.1.1.

- `tool_wrapper.ToolWrapper.run` crashed with `AttributeError:
  'CompletedProcess' object has no attribute 'pid'` after every
  subprocess call. `subprocess.run` returns a CompletedProcess after
  the child exits — the pid is gone, /proc/<pid> is gone, and
  per-child VmHWM is unrecoverable from this path. The mocked tests
  never exercised the real Popen→exit→stat path so the bug never
  fired. Set `peak_rss_mb=None` in this path with a NOTE referencing
  the Popen-based pipe path (Stage 3, where the read is still valid).

- `pileup_call._PILEUPCALLER_DATA_RE` used `[\d.]+` for the avg-coverage
  columns, which doesn't match scientific notation (e.g.,
  `3.6858006042296075e-2` emitted at low coverage). Switched to a full
  number regex `\d+(?:\.\d+)?(?:[eE][+-]?\d+)?`. Added a captured
  `pileupcaller_scientific.stderr` fixture from the smoke test as a
  regression test.

### Added (Day 7 — 2026-05-12)
- `pileup_aadr/coverage_impl.py` — Stage 0 diagnostic: wraps mosdepth
  with `--no-per-base --quantize 0:1:5:10:30 [--by REGIONS]` and parses
  `<prefix>.mosdepth.summary.txt` into `{per_chrom: {chrom: {length,
  bases, mean_coverage}}}`. Used to triage "BAM low-coverage globally"
  vs "BAM fine but AADR-1240k positions specifically lack coverage".
- `pileup_aadr/coverage_cmd.py` — click decorator: `BAM` positional,
  `--regions`/`--threads`/`--quantize`/`--json` options. Default TSV
  output; `--json` produces structured stdout for downstream tooling.
- `pileup_aadr/types.CoverageCliArgs` — frozen dataclass mirroring the
  click options 1:1 (4 fields + `bam`).
- `pileup_aadr/cli.py` — registers `coverage` on the root group.
  Removed the "Day 6-7 deferred" placeholder.
- `.github/workflows/ci.yml` — CI matrix: 6 cells (Python 3.11/3.12/
  3.13 × ubuntu-latest/macos-latest), `pip install -e '.[dev]'` then
  pytest. Separate `lint` job runs `ruff check` (strict) and `mypy`
  (informational — pandas-stubs gaps are not blocking until upstream
  stub coverage improves). Concurrency group cancels in-flight runs
  on superseded pushes.

### Tests added (Day 7; 208 total now)
- `test_coverage.py` — 6 tests: `_parse_mosdepth_summary` round-trips
  per-chrom dict; defensive skip on truncated rows; mocked-mosdepth
  TSV output round-trip; `--json` mode emits the same dict serialized;
  CLI root help lists `coverage`; `coverage --help` lists --regions/
  --threads/--json. Mosdepth itself is not required in the test env
  (ToolWrapper._resolve_binary is monkeypatched).

### Added (Day 6 — 2026-05-12)
- `pileup_aadr/types.py` — `ExtractCliArgs` frozen dataclass mirroring
  the `extract` subcommand's ~30 click options 1:1. Constructed via
  `ExtractCliArgs(**click_kwargs)` from the click decorator.
- `pileup_aadr/concurrency.py` — two context managers + helpers:
    - `output_lock(prefix)` — advisory `fcntl.flock` on
      `<prefix>.lock`. PID written to a SEPARATE `<prefix>.lock.holder`
      sidecar (H8 fix — releasing the lock no longer truncates the
      diagnostic the next contender reads). Contention raises
      `OutputLockHeldError` naming the holder PID.
    - `tempdir(*, base, keep_always, clean_on_crash)` — per-invocation
      tempdir with crash-survival semantics. Default: clean on success,
      retain on crash for forensics. Both directions overridable.
    - `warn_if_networked_fs(prefix)` — Linux-only stderr WARN when
      output is on NFS/SMB/CIFS where flock may be no-op'd.
- `pileup_aadr/output.py` — four writers + schema-version constants:
    - `write_pseudohaploid_sidecar` — `<prefix>.pseudohaploid.json`
      (consumed by pgen-samplebind); auto-injects schema_version.
    - `write_json_report` — schema-1 run summary with `tool` /
      `input` / per-stage `ExtractCounters` fields at top level (via
      `dataclasses.asdict`) / `output` / `config` blocks. No-lift
      fast path serializes Stages 1/2/4 as `null`.
    - `write_per_variant_tsv` — streaming 6-col TSV; constant memory
      regardless of variant count.
    - `write_stdout_summary` — human-readable multi-line block matching
      HLD §"Stdout summary" exactly: input/Stages 1-4/coverage report
      with 6-row autosome grid/output paths/wallclock.
- `pileup_aadr/extract_orch.py` — `run_extract(args)` orchestrator
  threading the 7-step HLD sequence: pre-flight (format detection +
  build detection + chain/ref resolution) → tool version probes →
  no-lift dispatch decision → `output_lock` + `tempdir` context →
  Stages 1-4 (or just 3 for fast path) → coverage gate evaluation →
  output writers → cleanup. H6 fix: output-prefix collision check
  happens INSIDE the locked region. C4 fix: `_run_stages` returns
  `(counters, rejoin_out)` so the populated sidecar dict + per-variant
  rows flow directly to `_write_outputs` (no rebuild from counters that
  loses `het_count`/`het_rate`).
- `pileup_aadr/extract_cmd.py` — replaced Day-1 NotImplementedError
  stub with a thin click wrapper that flips `--enable-baq` to the
  orchestrator's `no_baq=False`, constructs `ExtractCliArgs`, calls
  `run_extract`, and exits with the returned code.

### Tests added (Day 6; 202 total now)
- `test_concurrency.py` — 10 tests: lock acquire/release leaves 0-byte
  lock + removes holder sidecar; parent dir auto-created; real-
  subprocess contention raises with PID diagnostic; tempdir clean exit
  removes dir; default crash retention; `keep_always` retains on clean
  exit; `clean_on_crash` removes on crash; `pileup-aadr-` prefix; non-
  Linux fs detection returns None; local-FS path emits no NFS warning.
- `test_output.py` — 9 tests: sidecar JSON round-trip; schema_version
  injection; parent-dir creation; JSON report has tool block + counters
  fields at top level + correct values; no-lift JSON has null Stages
  1/2/4; per-variant TSV header + row count; iterator-based streaming;
  stdout summary lift-path renders all 4 stages; no-lift summary skips
  Stages 1/2/4 cleanly.
- `test_extract_orch.py` — 8 integration tests: end-to-end no-lift
  fast path writes all 4 artifacts; --report-json populates schema-1
  with null lift stages; --report-tsv populates 50-site streaming TSV;
  output-collision raises OutputExistsError; --overwrite replaces
  existing files; coverage gate failure raises CoverageGateFailure;
  warn-coverage threshold logs WARNING; full lift path dispatches to
  Stage 1 (mocked Picard, mocked pileup_call, mocked binary lookup).

### Added (Day 5 — 2026-05-12)
- `pileup_aadr/rejoin.py` — Stage 4 + no-lift fast-path finalizer.
  Two entry points share a `RejoinOutput` bundle (Stage4 counters,
  CoverageCounters, per-variant rows, sidecar dict):
    - `rejoin_aadr_frame()` — walks pileupCaller's per-variant `.geno`
      + `.snp` in lockstep, looks each rsID up in the AADR DataFrame,
      reads the `SwappedAlleles` flag from Picard's lifted-VCF INFO via
      pysam (`_build_swap_lookup`), applies `SWAP_DOSAGE` inversion
      (0->2, 1->1, 2->0, 9->9; involutive) when the flag fires, and
      streams the final EIGENSTRAT triplet to `<output_prefix>.{geno,
      snp,ind}` in AADR's hg19 frame. Defensive sanity check: lifted
      REF/ALT must match AADR REF/ALT modulo swap; mismatches drop with
      WARNING + counter. Per-chrom call counters + autosomal coverage
      gate computation alongside the main loop.
    - `_no_lift_fast_path_finalize()` — for AADR-build == BAM-build,
      pileupCaller's output IS the AADR-frame triplet. Copies `.geno`
      + `.snp` byte-for-byte via `shutil.copy2`, walks them once for
      counters, overrides `.ind` with user-supplied SEX (pileupCaller
      always writes SEX=U). `Stage4RejoinCounters.ref_alt_swap_count`
      and `allele_mismatch_drops` are 0 by definition.
- Module-level constants `GENO_HOM_REF` (0), `GENO_HET` (1),
  `GENO_HOM_ALT` (2), `GENO_MISSING` (9) — single ASCII chars matching
  pileupCaller's `.geno` on-disk encoding directly (no int lookups,
  no parse-and-reformat).
- PSEUDOHAPLOID sidecar JSON construction (consumed by pgen-samplebind):
  `pseudohaploid=1`, `het_count`, `non_missing_autosomal_count`,
  `het_rate`, `calling_mode="randomDiploid"`, `note` differentiating
  the lift vs no-lift paths.

### Tests added (Day 5; 175 total now)
- `test_rejoin.py` — 15 tests: SWAP_DOSAGE involutive correctness;
  no-swap passthrough writes EIGENSTRAT triplet; output `.snp` carries
  AADR's hg19 coords + Morgans (NOT pileupCaller's lifted hg38 coords);
  `.ind` 3-col TSV with user `--sex`; SwappedAlleles inverts dosage
  (0->2, 1->1, 2->0, 9->9 across all four geno chars in one run);
  swap-flag-but-alleles-don't-swap defensive drop + WARNING; no-swap-
  but-alleles-differ defensive drop + WARNING; rsID-in-pc-output-
  missing-from-AADR skipped + WARNING; malformed pileupCaller `.snp`
  row (!= 6 cols) → PileupAadrInternalError; coverage fraction =
  autosomal_calls / autosomal_aadr (sex chroms in per-chrom but not
  gated); `--report-tsv` per-variant rows populated with
  passthrough/swap/missing_call action labels; sidecar pseudohaploid
  classification with het_count + het_rate; no-lift fast-path triplet
  copy; no-lift sidecar note mentions "no-lift"; no-lift overrides
  pileupCaller's `.ind` with user SEX.

### Added (Day 4 — 2026-05-12)
- `pileup_aadr/transform.py` — Stage 2: `build_pileupcaller_snp_and_bed`
  reads Picard's lifted VCF and emits two artifacts:
    - **pileupCaller `.snp`** (6-col EIGENSOFT, AADR-numeric chrom in
      col 2 — chr1-22 → 1-22, chrX → 23, chrY → 24, chrM → 90; AADR
      rsID preserved from the AADR_RS INFO field with fallback to ID col).
    - **mpileup BED** (3-col 0-based, chr-prefixed to match modern hg38
      BAM @SQ headers).
  Alt-contig filter (default-on) drops alt/decoy contigs against the
  canonical-chrom regex `^(chr)?([0-9]{1,2}|X|Y|MT|M)$` derived from
  pileupCaller's parseSnpFile source — without it pileupCaller crashes
  mid-Stage-3 with an uncatchable Haskell SeqFormatException 5-10 minutes
  in. Defensive multi-allelic skip + numeric-chrom-map backstop covers
  the alt_contig_filter=False edge case.
- `pileup_aadr/pileup_call.py` — Stage 3: `run_pileup_call` builds the
  `samtools mpileup -B -q30 -Q30 -R -f <fasta> -l <bed> <bam>` →
  `pileupCaller --randomDiploid --seed N -f <snp> --sampleNames X
  --samplePopName Y -e <prefix>` pipe via `ToolWrapper.pipe`. Default
  threads cap at 4 (mpileup is BAM-seek-bound; verified empirically v2.1)
  with `no_thread_cap` opt-out. SIGPIPE handling: tolerates upstream
  exit 141 IFF downstream exited 0 (downstream is checked first); any
  other non-zero combination raises `ToolSubprocessError` with the
  stderr tail in the diagnostic. `parse_pileupcaller_stderr` extracts
  the structured 6-col TSV summary block (SampleName, TotalSites,
  NonMissingCalls, avgRawReads, avgDamageCleanedReads, avgSampledFrom);
  missing header or data line → `PileupAadrInternalError` with format-
  change diagnostic.

### Tests added (Day 4; 160 total now)
- `test_transform.py` — 9 tests: round-trip writes both files; numeric
  chrom encoding (1-22 / 23 / 24 / 90); chr-prefixed 0-based BED;
  AADR_RS INFO preferred over ID col; ID-col fallback when AADR_RS
  absent; alt-contig filter default drops chr*_random / chrUn_GL* /
  HLA-* contigs; alt_contig_filter=False with numeric-chrom-map
  backstop; empty-VCF empty-output path; output-dir auto-creation.
- `test_pileup_call.py` — 9 tests: stderr parser clean + missing-header
  + missing-data-line raises; mocked-pipe clean run populates
  Stage3CallCounters; downstream non-zero → ToolSubprocessError naming
  pileupCaller in `what`; upstream 141 + downstream 0 tolerated;
  upstream other non-zero → ToolSubprocessError naming samtools;
  thread-cap default applied + INFO-logged; `no_thread_cap=True`
  skips capping silently. ToolWrapper._resolve_binary monkeypatched
  to bypass on-PATH lookup (test env lacks samtools/pileupCaller).
- `tests/fixtures/stderr/pileupcaller_clean.stderr` — captured-from-v2.1
  pileupCaller 1.6.0.0 stderr summary block for parser unit tests.

### Added (Day 3 — 2026-05-12)
- `pileup_aadr/sites_vcf.py` — `build_sites_vcf` constructs a minimal
  VCF v4.2 from the AADR DataFrame, the shape Picard's LiftoverVcf
  expects (an EIGENSOFT `.snp` is not). Emits `##contig` lines
  matching the source build (hg19 by default; hg38 branch present
  for v0.2 forward-compat), `AADR_RS=<rsid>` INFO per row for
  Stage-4 lookup, and `chrN` form (normalizes AADR's numeric
  chrom). Palindrome filter (default-on) drops A/T + C/G ambiguous
  SNPs (~8% of biallelic 1240k); non-SNP filter (default-on) drops
  indels and non-ACGT alleles. Returns `Stage1InputFilters` counters.
  Custom text writer ~3x faster than pysam.VariantFile at 1.2M
  scale (verified at design time).
- `pileup_aadr/lift.py` (Stage 1 added on top of Day 2's chain
  resolution) — `lift_aadr_sites()` runs Picard LiftoverVcf via
  `ToolWrapper` with `RECOVER_SWAPPED_REF_ALT=true` +
  `WRITE_ORIGINAL_POSITION=true` + `WARN_ON_MISSING_CONTIG=true`
  + `MAX_RECORDS_IN_RAM=100000`. JVM heap (`-Xmx<picard_mem>`,
  default 3g) routed via H9-fixed `jvm_args` parameter (not
  positional args). `parse_picard_stderr()` extracts 5 counters
  via REQUIRED + OPTIONAL regex pattern split (M14 fix — `swapped`
  is OPTIONAL, defaults 0 if Picard ever omits it). REQUIRED-pattern
  miss → `PileupAadrInternalError` with reproducible "Picard format
  changed" diagnostic. `parse_rejected_vcf()` categorizes by
  FILTER (NoTarget / MismatchedRefAllele /
  IndelStraddlesMultipleIntervals / SwappedAlleles / other);
  empty path → all-zero (defensive). Yield gate: fail
  (`LiftoverYieldError`) below 70% (default), warn (stderr WARN +
  JSON flag) below 95% (default).

### Tests added (Day 3; 142 total now)
- `test_sites_vcf.py` — 13 tests: round-trip from Day-1 chr22
  fixture (50 sites, 1 palindrome → 49 written), VCF v4.2 header
  correctness (##contig length matches HG19_CHROM_LENGTHS, AADR_RS
  INFO present, 8-col CHROM line), chrN normalization,
  palindrome filter on/off, non-SNP filter (indels + N + IUPAC
  codes), unrecognized-chrom DEBUG-skip path, sort order
  (CHROM_ORDER then pos), aadr_build switch (hg19 vs hg38
  ##contig table), empty-input header-only path.
- `test_lift_stage1.py` — 9 tests: parse_picard_stderr against
  captured Picard 3.3.0 fixtures from v2.1 verification (clean
  100-site run + partial-yield 5000-site run), optional-swap
  default, missing-required-pattern raises; parse_rejected_vcf
  empty + categorization paths; lift_aadr_sites mocked-Picard
  integration covering clean run (100% yield), low-yield raise
  (50% < 70% gate), warn-only path (80% < 95% warn but ≥ 70% fail).
- `tests/fixtures/stderr/picard_clean.stderr` and
  `picard_partial_yield.stderr` — captured-from-v2.1 Picard 3.3.0
  output for stderr-parser unit tests.

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
