# DEVELOPMENT.md

A guide for someone picking up the baton on `pileup-aadr`. Read this in order
the first time. After that, jump to the section relevant to your change.

For user-facing install + usage, see [README.md](README.md). For the
day-to-day dev loop (test invocation, lint, release), see
[CONTRIBUTING.md](CONTRIBUTING.md). This file is about *understanding the
code well enough to extend it*.

> **Note on current state.** The repo is at **v0.5.0** (published to PyPI).
> This file is the source of truth on architecture. Exact line-number
> references below may drift between releases — trust the described
> contracts over the line numbers.

---

## 1. What the tool does (in 90 seconds)

`pileup-aadr extract <BAM> <AADR.snp> -o <prefix>` produces an EIGENSTRAT
triplet (`.geno`/`.snp`/`.ind`) of pseudohaploid genotypes for a user's
WGS BAM at every site in an AADR `.snp` panel. Downstream, the triplet is
sample-bound to the AADR cohort by `pgen-samplebind` and consumed by
AdmixTools 2 (`qpAdm`, `qpgraph`).

It replaces the 5-step manual pipeline that every ancient-DNA personal-WGS
study currently reimplements: lift AADR sites → sites VCF → pileupCaller
`.snp` → mpileup + `pileupCaller --randomHaploid` → rejoin to AADR's hg19
frame by rsID. Each step has historically-tribal failure modes; wrapping
all five in one tool with proper instrumentation removes the most
failure-prone manual link in any aDNA pipeline.

The work happens in four stages — **lift, transform, call, rejoin** —
gated by a coverage check at the end. When the BAM and AADR are already
in the same build (e.g. a hg19 BAM + hg19 AADR), Stages 1/2/4 are skipped
entirely (the **no-lift fast path**); pileupCaller's output IS the final
triplet. Otherwise Picard lifts AADR hg19 → hg38 with
`RECOVER_SWAPPED_REF_ALT=true`, transform builds pileupCaller's input
files, mpileup-pipe-pileupCaller emits hg38-frame genotypes, and rejoin
walks the result row-by-row, using rsID to land each call back in AADR's
hg19 frame, inverting dosage where Picard set the `SwappedAlleles` INFO
flag.

Three companion subcommands ship in the same CLI: `validate` (a 16-result
pre-flight that's cheap to run before a 30-min `extract`), `coverage` (a
mosdepth wrapper for triaging coverage-gate failures), and `inspect` (a
pure-Python AADR `.snp` summary).

---

## 2. The 4-stage pipeline (and the no-lift fast path)

```
AADR .snp (hg19)              user BAM (hg19 or hg38)
       │                            │
       │ Stage 1: lift.py
       │   sites_vcf.build_sites_vcf  → minimal VCF v4.2 (palindrome + non-SNP filtered)
       │   lift.lift_aadr_sites       → Picard LiftoverVcf RECOVER_SWAPPED_REF_ALT=true
       │   (skipped on no-lift fast path)
       ▼
   lifted VCF + SwappedAlleles INFO flags
       │
       │ Stage 2: transform.py
       │   build_pileupcaller_snp_and_bed → .snp (numeric chrom) + .bed (chr-prefixed)
       │   alt-contig filter (chr*_*_alt) drops Picard's strand-jumped hits
       │   (skipped on no-lift fast path)
       ▼
   .snp + BED in target build
       │
       │ Stage 3: pileup_call.py
       │   samtools mpileup -B -q 30 -Q 30 -R -f <fa> -l <bed>  |
       │     pileupCaller --randomHaploid   (default; --calling-mode selects)
       │   (always runs; ~25–40 min on a 33× WGS at 1240k)
       ▼
   pileupCaller EIGENSTRAT triplet (target build)
       │
       │ Stage 4: rejoin.py
       │   rejoin_aadr_frame → per-row walk, AADR-frame keys, dosage-flip
       │     where SwappedAlleles, streams final .geno/.snp/.ind
       │   _no_lift_fast_path_finalize (fast path) → trivial pass-through
       ▼
   final EIGENSTRAT triplet + <prefix>.pseudohaploid.json sidecar
       │
       │ Coverage gate (extract_orch._evaluate_coverage_gate)
       │   absolute autosomal call count vs --min-coverage / --warn-coverage
       │   panel-aware: skipped with N/A for chrY-only / chrM-only / sex-only
       ▼
   exit 0 (or exit 1 on FAIL gate)
```

**The no-lift fast path is the simplifier.** When `--bam-build == --aadr-build`,
the orchestrator branches to a much shorter sequence: Stage 3 only,
followed by a copy-through finalize. Stages 1/2/4 serialize as `null` in
the JSON report.

---

## 3. Codebase tour

The package is **~6,790 LOC across 28 modules** under `pileup_aadr/`.
Names to know on day one:

| Module | Job |
|---|---|
| `cli.py` | Click root group; outermost `try/except PileupAadrError` formatter; maps to exit codes; calls `configure_logging` once |
| `extract_cmd.py` | Click decorator + **31** options for `extract`; constructs `ExtractCliArgs`, delegates to orchestrator |
| `extract_orch.py` | **Orchestrator.** Pre-flight → tool versions → lock+tempdir → stages → gates → writers |
| `validate_cmd.py` / `validate_impl.py` | `validate` subcommand (16 results across pre-flight checks) |
| `coverage_cmd.py` / `coverage_impl.py` | `coverage` subcommand (mosdepth wrapper) |
| `inspect_cmd.py` / `inspect_impl.py` | `inspect` subcommand (pure-Python AADR `.snp` summary) |
| `format_detect.py` | BAM/CRAM detection, build detection (chr1 + chr20 anchors), sample-name resolution, AADR parser, `classify_aadr_chrom_set`, `BuildOverride = Literal["hg19", "hg38", "auto"]` |
| `sites_vcf.py` | Pre-Stage-1: AADR DataFrame → minimal v4.2 sites-VCF with palindrome + non-SNP filters |
| `lift.py` | **Stage 1** + chain resolution (3-tier: CLI → env → bundled) + bundled chain SHA verification + Picard stderr/rejected-VCF parsers |
| `transform.py` | **Stage 2** |
| `pileup_call.py` | **Stage 3** + pileupCaller stderr summary parser |
| `rejoin.py` | **Stage 4** + `SWAP_DOSAGE` table + streaming `.geno`/`.snp`/`.ind` writes + `_no_lift_fast_path_finalize` |
| `tool_wrapper.py` | Subprocess plumbing — `ToolWrapper`/`ToolSpec`/`ToolRunResult`; one `*_SPEC` per binary (samtools/pileupCaller/picard/java/mosdepth) |
| `ref_resolve.py` | Target FASTA resolution (CLI/env/BAM @PG) + chr1-length build verification |
| `dict_resolve.py` | `.dict` find-or-auto-generate (sibling-of-FASTA → user cache fallback) |
| `concurrency.py` | `output_lock(prefix)` flock + `tempdir(...)` lifecycle + networked-FS warning |
| `counters.py` | 7 frozen per-stage counter dataclasses + the mutable `ExtractCounters` aggregate |
| `output.py` | Sidecar JSON, run-level JSON report, streaming per-variant TSV, stdout summary writers |
| `errors.py` | Exception hierarchy keyed by exit code (1/2/3/4) with `__init_subclass__` enforcement |
| `types.py` | Frozen `ExtractCliArgs` / `CoverageCliArgs` (CLI kwargs view, decoupled from Click) |
| `chrom_lengths.py` | Hardcoded hg19/hg38 chrom-length tables |
| `logging_config.py` | stderr logger; `$PILEUP_AADR_JSON_LOGS=1` switches to JSON Lines |
| `data/` | Bundled `hg19ToHg38.over.chain.gz` (~223 KB) + `.sha256` sidecar |

> **Naming gotcha.** The four stages live in semantically named files
> (`lift.py`, `transform.py`, `pileup_call.py`, `rejoin.py`), NOT in
> `extract_stage1.py`/`stage2.py`/etc.

---

## 4. Walking through `pileup-aadr extract`

```mermaid
sequenceDiagram
    autonumber
    participant U as User shell
    participant CLI as cli.py
    participant CMD as extract_cmd.py
    participant ORCH as extract_orch<br/>run_extract
    participant STAGES as Stage modules<br/>(lift/transform/<br/>pileup_call/rejoin)
    participant OUT as output.py

    U->>CLI: pileup-aadr extract BAM AADR -o PREFIX
    CLI->>CLI: main() — configure_logging;<br/>outer try/except PileupAadrError
    CLI->>CMD: dispatch to extract() Click cmd
    CMD->>CMD: collect 31 options; invert<br/>--enable-baq → no_baq;<br/>build frozen ExtractCliArgs
    CMD->>ORCH: run_extract(args)
    ORCH->>ORCH: pre-flight: format, build, sample-name,<br/>AADR parse, chain, FASTA, .dict
    ORCH->>ORCH: probe binary versions<br/>(ToolWrapper.version per SPEC)
    ORCH->>ORCH: acquire output_lock(prefix);<br/>create tempdir
    alt Full path: BAM build ≠ AADR build
        ORCH->>STAGES: Stage 1 — sites_vcf + Picard LiftoverVcf
        ORCH->>STAGES: Stage 2 — pileupCaller .snp + BED
        ORCH->>STAGES: Stage 3 — samtools mpileup → pileupCaller
        ORCH->>STAGES: Stage 4 — rejoin_aadr_frame (streaming<br/>writes the final triplet directly)
    else No-lift fast path: BAM build == AADR build
        ORCH->>STAGES: Stage 3 only
        ORCH->>STAGES: _no_lift_fast_path_finalize<br/>(copy-through)
    end
    ORCH->>ORCH: _evaluate_coverage_gate;<br/>_finalize_counters
    ORCH->>OUT: write_pseudohaploid_sidecar (always)
    ORCH->>OUT: write_run_report_json (if --report-json)
    ORCH-->>CMD: exit_code (int)
    CMD-->>CLI: ctx.exit(exit_code)
    CLI-->>U: stable exit code 0/1/2/3/4
```

**Ownership rules:**
- The orchestrator owns the lock + tempdir lifecycle, JSON report assembly,
  gate evaluation, and the no-lift branch.
- Stage modules own their stage's work and their stage's `Counters`
  dataclass. They never touch the lock, never write to the output prefix.
  **Stage 4 is the only exception** — it streams the final triplet
  directly because the alternative is materializing 1.2M rows in RAM.

**The other three subcommands** are much shorter:
- `validate` (`validate_impl.run_validate`) → runs 16 checks, emits TSV
  or JSON, returns 0/1.
- `coverage` (`coverage_impl.run_coverage`) → wraps mosdepth, parses its
  per-chrom output, emits TSV or JSON.
- `inspect` (`inspect_impl.run_inspect`) → pure-Python AADR `.snp` parse +
  summary (panel size, allele distribution, palindrome %, `chrom_set`).

---

## 5. Key contracts (don't break these)

### 5.1 Exit codes (stable across versions)

| Code | Bucket | When |
|---|---|---|
| 0 | Success | Possibly with `WARN`-level gate flags |
| 1 | Soft validation | `LiftoverYieldError`, `CoverageGateFailure` — data-quality issues |
| 2 | I/O | Chain/FASTA/BAM not found, subprocess crash, lock held, `.bai` missing, etc. |
| 3 | Invariant | Build mismatch, AADR malformed, defensive sanity, uncaught exception |
| 4 | Usage / dependency | Bad CLI args, missing or wrong-version external binary |

Workflow managers code against these codes. Don't add a new bucket.
If you add a new error class, it inherits from `PileupAadrError` and
declares `exit_code` ∈ {1, 2, 3, 4} — `__init_subclass__` enforces it.

### 5.2 Error class taxonomy (`pileup_aadr/errors.py`)

**25 named subclasses** of `PileupAadrError`, ordered in the file BY
EXIT CODE (not alphabetical — `pyproject.toml` `RUF022` ignore for
readability of the taxonomy). Note: `ToolNotFoundError` has 5 concrete
subclasses (one per binary). Every exception ships with a `What:` /
`Why:` / `Fix:` triple in its constructor; `errors.format_error`
renders them. Renaming classes (e.g. adding an `Error` suffix to
`BAMSampleNameAmbiguous`) breaks the spec contract — the
`pyproject.toml` `N818` ignore exists for this reason.

### 5.3 JSON report schema (v2, additive-only)

Top-level keys, built by `output.write_run_report_json`:

```
schema_version
tool { name, version, tool_versions{} }
input { bam_path, bam_format, bam_build, bam_sample_name_source,
        aadr_snp_path, aadr_build, aadr_input_rows,
        chain_path, chain_sha256, target_fasta_path }
stage_1_lift     { ... } | null   # null on no-lift fast path
stage_2_transform { ... } | null  # null on no-lift fast path
stage_3_call     { wallclock_seconds, pileupcaller_summary { ... } }
stage_4_rejoin   { ... } | null   # null on no-lift fast path
coverage         { non_missing_autosomal_calls, coverage_fraction,
                   coverage_warning, per_chrom_call_count{} }
gates            { liftover_yield: PASS|WARN|FAIL|N/A, coverage: same }
output           { prefix, geno_bytes, snp_bytes, ind_bytes,
                   pseudohaploid_sidecar }
config           { resolved settings }
wallclock_total_seconds
exit_code
```

**Additive-only across versions.** Two consumers parse this:
- ancestry-pipeline-tool's gate node reads `coverage.non_missing_autosomal_calls`
  and `gates.coverage`
- pgen-samplebind reads the separate `<prefix>.pseudohaploid.json` sidecar
  for `pseudohaploid: 1`, `het_count`, `het_rate`, `non_missing_autosomal_count`

Renaming or removing a field is a breaking change. Adding a field
underneath an existing key is fine; permissive consumers tolerate
unknowns.

### 5.4 PSEUDOHAPLOID sidecar (`<prefix>.pseudohaploid.json`)

Separate file from the run report. Read by `pgen-samplebind` at merge
time to skip per-sample heterozygosity recomputation. Required fields:
`pseudohaploid` (1 for the pseudo-haploid modes `randomHaploid`/`majorityCall`,
0 for the diploid `randomDiploid` escape hatch), `het_count`,
`non_missing_autosomal_count`, `het_rate`, `source: pileup-aadr-extract`,
`calling_mode` (the `--calling-mode` value, default `randomHaploid`).

### 5.5 External-binary version contract

Probed at orchestrator startup via `tool_wrapper.ToolWrapper(SPEC).version()`
before any stage runs. Cached into `tool.tool_versions` of the JSON
report. Failure modes: missing → `ToolNotFoundError` (exit 4), wrong
version → `ToolVersionError` (exit 4), unparseable output →
`PileupAadrInternalError` (exit 3, defensive sanity). Min versions:
samtools ≥ 1.16, pileupCaller ≥ 1.6.0, picard ≥ 3.0, java ≥ 11,
mosdepth ≥ 0.3.6.

---

## 6. Architecture details (the internals you'll touch)

### 6.1 `tool_wrapper.py` — the subprocess discipline

Every external binary call goes through a `ToolWrapper(spec)` instance.
The 5 `*_SPEC` constants (`SAMTOOLS_SPEC`, `PILEUPCALLER_SPEC`,
`PICARD_SPEC`, `JAVA_SPEC`, `MOSDEPTH_SPEC`) carry: binary name,
`version_args` (e.g., `["--version"]`), version regex, min version.

- `version()` — invokes the version probe. **Uses `subprocess.run`
  directly, NOT `self.run`**, because (1) the version probe predates
  the tempdir/stderr-disk plumbing and (2) output is bounded < 1 KB.
- `run(args, capture_stderr_to: Path, ...)` — runs a tool with stderr
  always written to `capture_stderr_to` (tail-readable for diagnostic).
  Returns a `ToolRunResult` carrying `stdout`, `stderr_path`,
  `stderr_text` (populated only when explicitly read), elapsed time,
  RSS measurement.
- `pipe(left_args, right_args, ...)` — the Stage 3 idiom: builds a
  `samtools mpileup … | pileupCaller …` Popen pair with SIGPIPE
  handling.

The `stderr_path` discipline is load-bearing: invariant #7 (no
unbounded stderr) depends on the orchestrator passing a tempdir path
into every `run()` call.

### 6.2 Resolution chains — chain file + reference FASTA

Both follow the same 3-tier pattern. The orchestrator resolves once
during pre-flight and threads the resolved `Path` into stage modules.

**Chain (`lift.chain_file_path`):**
1. `--chain PATH` (explicit CLI) — used as-is; SHA verification only if
   `--strict-chain-sha`.
2. `$PILEUP_AADR_CHAIN_DIR/hg19ToHg38.over.chain.gz` — same.
3. `lift.get_bundled_chain_path()` — `importlib.resources.files("pileup_aadr") / "data" / "hg19ToHg38.over.chain.gz"`. **SHA always verified** against the bundled `.sha256` sidecar; mismatch → `ChainFileSHAError` with reinstall guidance.

`--insecure-chain` skips SHA verification entirely with a stderr WARNING.

**Target FASTA (`ref_resolve.resolve_ref_fasta`):**
1. `--ref-fasta PATH` (explicit CLI).
2. `$PILEUP_AADR_REF_DIR/<build>.fa` (e.g., `hg38.fa`).
3. BAM `@PG` autodetect — pysam reads `bam.header.to_dict().get('PG', [])`
   for path hints embedded by upstream aligners.
4. Otherwise → `ReferenceFastaNotFound` (exit 2).

`ref_resolve.verify_fasta_matches_bam_build` runs after resolution:
checks the FASTA's chr1 length matches the BAM's build before Picard
spends 30+ seconds rejecting most sites with `MismatchedRefAllele`.

**FASTA `.dict` (`dict_resolve.ensure_target_fasta_dict`):** sibling-of-FASTA
first, then user cache (`~/.cache/pileup-aadr/dicts/<basename>.dict`),
then auto-generate via `picard CreateSequenceDictionary` (~23s for
hg38, cached).

### 6.3 Output lock + tempdir (`concurrency.py`)

- `output_lock(prefix)` — opens `<prefix>.lock` with
  `fcntl.LOCK_EX | LOCK_NB`. Lock release is via fd close on context
  exit (file left in place — avoids race with unlink). A separate
  `<prefix>.lock.holder` PID sidecar is written on acquire, removed on
  release; reading it can race but is best-effort diagnostic.
- `tempdir(base, keep_always, clean_on_crash)` — creates
  `pileup-aadr-XXXXXX/` under `$TMPDIR` (or `--tempdir`). Default:
  clean on success, RETAIN on crash with path logged at ERROR.
  `--keep-tempdir` retains always; `--clean-tempdir-on-crash` cleans
  even on crash (CI/container use).
- `warn_if_networked_fs()` parses `/proc/mounts` (Linux only) and
  emits a stderr WARNING for nfs/cifs mounts (where flock may be a
  silent no-op).

The orchestrator entry sequence is: `output_lock` first, then
`tempdir` inside. The output-prefix existence check happens **inside**
the lock (closes a race window where two concurrent invocations could
both see "no existing output" and both proceed).

### 6.4 Stage-4 streaming (`rejoin.rejoin_aadr_frame`)

Why streaming matters: a 1.2M-site EIGENSTRAT matrix is small but
non-trivial. The design constraint is to keep memory flat regardless
of panel size.

How it works:
1. **Pre-build the swap lookup** — walk the lifted VCF once, build a
   `{rsid: bool}` dict recording which rsIDs Picard flagged with
   `SwappedAlleles`.
2. **Zip pileupCaller's `.geno` + `.snp` outputs** — open both,
   iterate as `for geno_line, snp_line in zip(pc_geno, pc_snp):`.
3. **Per row** — read the rsID from the snp line; look up the AADR row
   via `aadr_df.loc[rsid]`; if `swap_lookup.get(rsid, False)`, apply
   `SWAP_DOSAGE[geno_char]` (the table at the top of `rejoin.py`:
   `{"0": "2", "1": "1", "2": "0", "9": "9"}`); write the AADR-frame
   row immediately to the final `.geno`/`.snp`.

The streaming write IS the rejoin — splitting writer-from-loop would
require iterators with row context. `--report-tsv` per-variant rows
accumulate only when `emit_per_variant_rows=True`.

### 6.5 `validate` — the 16 emitted results

`validate_impl.run_validate` builds a `list[CheckResult]` in this
order. Each row emits one result; the function returns exit 1 if any
result is `FAIL`, else exit 0.

| # | Check | Source |
|---|-------|--------|
| 1 | AADR `.snp` parseable + no dup rsIDs | `_aadr_parse_pass` (inline) |
| 2 | AADR build detectable | `detect_aadr_build` (inline) |
| 3 | BAM index present + readable | `_bam_index_check` |
| 4 | BAM build detectable | `detect_bam_build` (inline) |
| 5 | samtools version | `_check_tool(SAMTOOLS_SPEC)` |
| 6 | picard version (or SKIP on no-lift) | `_check_tool(PICARD_SPEC)` |
| 7 | java version (or SKIP on no-lift) | `_check_tool(JAVA_SPEC)` |
| 8 | pileupCaller version | `_check_tool(PILEUPCALLER_SPEC)` |
| 9 | samtools mpileup flag probe | `_check_tool_flags_samtools_mpileup` |
| 10 | pileupCaller flag probe | `_check_tool_flags_pileupcaller` |
| 11 | mosdepth flag probe | `_check_tool_flags_mosdepth` |
| 12 | picard LiftoverVcf flag probe (or SKIP on no-lift) | `_check_tool_flags_picard_liftover` |
| 13 | Chain file resolves + SHA-verified | `_check_chain` |
| 14 | Target FASTA findable + build matches | `_check_ref_fasta` |
| 15 | FASTA `.dict` present or generatable (or SKIP on no-lift) | `_check_target_fasta_dict` |
| 16 | Output prefix not lock-held + parent writable (or SKIP if no `-o`) | `_check_output_prefix` |

> The `run_validate` docstring still says "10 pre-flight checks" —
> that's the original spec count. #9–12 (flag probes) were added in
> v0.2 after [issue #2](https://github.com/carstenerickson/pileup-aadr/issues/2)
> burned a customer with a wrong-flag bug, and #15 was split out from
> the original FASTA check.

### 6.6 Logging

Set up exactly once by `cli.main` → `configure_logging(level=…)`.
Pattern: every module declares `log = logging.getLogger(__name__)`.

Default format (human-readable):
```
2026-05-13 10:23:41 [INFO] pileup_aadr.lift: Wrote lifted VCF: /tmp/.../lifted.vcf
```

JSON Lines (enable via `PILEUP_AADR_JSON_LOGS=1`):
```json
{"ts": "2026-05-13T17:23:41.123+00:00", "level": "INFO", "logger": "pileup_aadr.lift", "msg": "Wrote lifted VCF: /tmp/.../lifted.vcf"}
```

The JSON variant exists for ancestry-pipeline-tool's stderr-streams-to-disk
discipline (parses without regex). Extra fields attached via
`logger.info(..., extra={…})` flow through. Ruff `T20` forbids `print()`
in module code; only the orchestrator may write stdout (via
`sys.stdout.write` for the human summary block).

---

## 7. Hard invariants (the "don't break" list)

### Correctness (data integrity)

1. **`RECOVER_SWAPPED_REF_ALT=true`** drives Stage 1 lift; the
   `SwappedAlleles` INFO flag is the input to Stage 4's dosage flip.
   Without it, ~3% of sites silently get the wrong dosage.
2. **The default `--randomHaploid` (and `--majorityCall`) output is
   pseudohaploid by construction** (0% het, matches the pseudo-haploid
   AADR panel). The sidecar JSON records `calling_mode` + `pseudohaploid`;
   there's no 4th `.ind` column for it because convertf/mergeit silently
   drop unknown columns. The `--randomDiploid` escape hatch produces
   diploid (het-bearing) calls and is recorded as `pseudohaploid=0` so the
   downstream f2 consumer never treats diploid data as pseudo-haploid.
3. **Duplicate-rsID rejection at startup** — `format_detect.parse_aadr_snp`
   raises `AADRDuplicateRsidError` (exit 3) on dupes.
4. **`--sample-name` ALWAYS overrides BAM `@RG SM:`** — even when they
   disagree; emits a stderr INFO line for audit.
5. **`Stage1LiftCounters.swapped_alleles_count >= Stage4RejoinCounters.ref_alt_swap_count`**
   — Stage 1 is an upper bound on swaps; Stage 4's defensive sanity
   check drops sites where lifted REF/ALT don't actually swap to AADR's
   REF/ALT (those count toward `Stage4.allele_mismatch_drops`).
   Asserted in `tests/integration/test_real_binaries.py:124-130`.

### I/O discipline

6. **Bundled chain SHA verification at startup** — always-on (~2 ms
   cost). Mismatch → `ChainFileSHAError` with reinstall guidance.
7. **Subprocess stderr is always written to a per-tool log file**
   (`ToolRunResult.stderr_path`). The orchestrator routes these into
   the tempdir. Never collect stderr into unbounded memory.
8. **Streaming Stage-4 writes** — `rejoin_aadr_frame` zips pileupCaller
   `.geno` + `.snp` and writes the AADR-frame triplet inside the same
   loop. No materialization of the 1.2M-row matrix in RAM. See §6.4
   for the architecture.

### Concurrency

9. **Output-prefix existence check happens INSIDE the lock** to close
   a race with concurrent runs.
10. **`_check_output_prefix` in `validate` does NOT acquire the lock**
    — would give false confidence to a concurrent extractor. It only
    inspects the `<prefix>.lock.holder` PID sidecar.

### Compatibility & style

11. **Stable exit codes.** Workflow managers depend on the 0/1/2/3/4
    contract. Don't add a 5.
12. **No `print()` in module code.** Only the orchestrator may write
    stdout (via `sys.stdout.write` for the summary block); everything
    else logs through the stdlib `logging` module. Ruff `T20` enforces.

---

## 8. Gotchas (things that have bitten before)

### CLI quirks

- **The BAQ flag is doubly inverted.** `--enable-baq` (CLI flag) →
  `no_baq=True` (kwarg) → omits `-B` from `samtools mpileup` (which is
  the *disable* BAQ flag). The CLI default is BAQ-disabled (matches
  pileupCaller's recommended cmdline). The mapping is at
  `extract_cmd.py:233-242`; the full-chain f2 invariant test caught a
  sign-error here in v0.2.
- **`--quiet` and `--verbose` are root-group options.**
  `pileup-aadr extract --quiet …` rejects with "no such option"; use
  `pileup-aadr --quiet extract …`.
- **`samtools mpileup` has no `-@` flag, ever.** v0.1.0–0.1.1 passed it
  erroneously and burned [issue #2](https://github.com/carstenerickson/pileup-aadr/issues/2).
  Current code keeps the `--threads` CLI flag for back-compat but only
  ever passes it to pileupCaller; samtools is invoked single-threaded.
  A WARN log fires if `--threads > 1` so users know it's a no-op for
  mpileup.

### Naming & style

- **Error class names don't end in `Error`** (some do, some don't).
  Names like `BAMSampleNameAmbiguous` and `CoverageGateFailure` are
  pinned by the spec; ruff `N818` is disabled for this. Don't rename.
- **`errors.py` `__all__` ordering is by exit code, not alphabetical.**
  Ruff `RUF022` is disabled for this file specifically.
- **Benchmarks use `bench_*.py`, not `test_*.py`.** They live under
  `benchmarks/` (excluded from `testpaths`); `python_files` config
  matches both prefixes so explicit `pytest benchmarks/` collects them.
- **mypy `strict = true` is enforced in CI** (the lint job fails on any
  error — it is no longer informational). The tree is effectively
  `# type: ignore`-free: pandas `itertuples()` scalar unions are handled
  with `cast`, Click command callbacks make `ctx` positional-only, and
  pysam header reads go through `.to_dict()`. `pandas-stubs` is **pinned**
  in the dev extras (its `itertuples`/column typing shifts between
  releases, so an unpinned floor makes `mypy` non-reproducible between
  local and CI); `types-psutil` is declared for the same reason.

### Streaming & data shape

- **`ExtractCounters` is mutable but leaf counters are frozen.** The
  orchestrator replaces `coverage` and `gates` via `dataclasses.replace`
  after gate evaluation, not field-mutation.
- **Stage 4 writes the final triplet** — see §6.4.
- **chr1 lengths between hg19 (249,250,621) and hg38 (248,956,422)
  differ by only 294 KB.** Build detection uses *closest-match*, not
  first-match-wins; a ±1 Mb tolerance would tag every hg38 BAM as hg19.
  This was [issue #1](https://github.com/carstenerickson/pileup-aadr/issues/1).
  chr20 is the fallback anchor (1.4 Mb gap, more forgiving).

### pysam quirks

- **`pysam.AlignmentFile` open mode depends on what you're doing.**
  `format_detect.detect_bam_format` uses `"rb"` because it's
  specifically testing whether the file IS a BAM (the magic-byte path).
  All other call sites (`detect_bam_build`, `detect_bam_sample_name`)
  use `"r"` because that auto-detects BAM vs CRAM. Don't `"rb"` in the
  general path — v0.1 did and crashed on CRAM.
- **`@RG SM:` extraction** raises `BAMSampleNameAmbiguous` only if
  multiple *distinct* values exist. A single SM repeated is fine;
  the filename-stem fallback is silent.

---

## 9. Common tasks (where to start)

### Add a new CLI flag to `extract`
1. Add the `@click.option(...)` decorator in `extract_cmd.py:23-213`
   (the option block).
2. Add the matching field to `types.ExtractCliArgs` (frozen dataclass).
3. If the flag affects a stage, thread it through `extract_orch.run_extract`
   to that stage's call site.
4. Stage modules consume from their function signature, not a global.
5. Add a unit test in `tests/test_extract_cmd.py` covering CLI default
   + flag-set behaviour.
6. Add a CHANGELOG entry (`## [Unreleased]` heading or a new
   `## [vX.Y.Z]` heading if you're cutting a release).
7. Add to README's "Notable extract options" list if user-facing.

### Add a new error class
1. Pick the exit code (1/2/3/4).
2. Add the class to `pileup_aadr/errors.py` in the right exit-code group
   (file is ordered by exit code intentionally; don't sort).
3. `__init_subclass__` will reject anything missing `exit_code`.
4. Use `What:` / `Why:` / `Fix:` in the constructor-built message;
   `errors.format_error` renders them.
5. Add to `__all__` (also exit-code-ordered).
6. Add a test in `tests/test_errors.py` covering the message contract —
   every named exception's `format_error` output must contain
   `What:`/`Why:`/`Fix:`.
7. Add a CHANGELOG entry.

### Add a new `validate` check
1. Write `_check_<name>(args, ...) -> CheckResult` in `validate_impl.py`.
2. Wire it into `run_validate`'s sequence (currently 16 emitted results;
   see §6.5 for the ordered list).
3. Decide PASS/WARN/FAIL/SKIP semantics; FAIL contributes to exit 1.
4. Add a unit test in `tests/test_validate_cmd.py`.
5. Update the count in §6.5 of this file.
6. Add a CHANGELOG entry.

### Add a new field to the JSON report
1. Add the field to the relevant counter dataclass in `counters.py`.
2. The serializer in `output.py` uses `dataclasses.asdict()` — no manual
   plumbing required IF the field belongs in an existing stage's counter.
3. If it belongs at top-level (not inside a stage), update
   `output.write_run_report_json` directly.
4. Add to `tests/test_output.py` covering presence + type.
5. **Additive only** — don't rename existing fields. Permissive
   consumers (ancestry-pipeline-tool, pgen-samplebind) tolerate unknowns
   but break on missing-known-fields.
6. Add a CHANGELOG entry.

### Bump the bundled chain file
1. Replace `pileup_aadr/data/hg19ToHg38.over.chain.gz`.
2. Recompute SHA256 → write to `.sha256` sidecar.
3. The startup verifier picks it up automatically (~2 ms cost).
4. Run `pytest tests/test_lift_chain.py` and the full-chain f2
   invariant test (`tests/integration/test_lld19_full_chain.py`, Layer B
   with a real BAM) before tagging.
5. Add a CHANGELOG entry noting the chain source + date.

### Touch Stage 4 (rsID rejoin or SwappedAlleles handling)
1. **This is the most error-prone link in the chain.** Read §6.4 first.
2. The `SWAP_DOSAGE` table at the top of `rejoin.py` is the canonical
   answer for inverting EIGENSTRAT dosage characters.
3. The streaming write IS the rejoin — see §6.4 for the architecture.
4. The full-chain f2 invariant test defends correctness end-to-end:
   `tests/integration/test_lld19_full_chain.py` runs with `max_dev < 1e-9`
   against a frozen baseline. Re-run Layer B (env-gated on a real BAM)
   after any Stage 4 change.

---

## 10. Test discipline

- **Unit tests** live in `tests/test_<module>.py`, mock all subprocess
  calls via `pytest-mock`. CI matrix is 6 cells
  (Python 3.11/3.12/3.13 × Linux/macOS).
- **Integration tests** live in `tests/integration/`. They invoke real
  binaries and skip cleanly when toolchain absent. CI's `bio-tools`
  job sets up the toolchain via bioconda; locally see CONTRIBUTING.md.
- **Total: 313 tests collected** by plain `pytest`. The unit subset
  runs in ~5s; integration tests auto-skip without their toolchain, so
  the total runtime with no binaries on PATH is also ~5s.
- **The full-chain f2 invariant** (`tests/integration/test_lld19_full_chain.py`)
  is the integration test you must NOT break: it asserts the full chain
  pileup-aadr → pgen-samplebind → AT2 `extract_f2` produces f2 numbers
  stable to `max_dev < 1e-9` of a frozen baseline. Layer A runs in
  ~10s with bundled fixtures; Layer B requires a real BAM via
  `PILEUP_AADR_LLD19_BAM`.
- **Markers** (declared in `pyproject.toml`): `slow` (>10s; skipped by
  default), `requires_picard`, `requires_samtools`, `requires_pileupcaller`,
  `requires_mosdepth`. `addopts = "--strict-markers"` rejects unknowns.
- **Benchmarks** (`benchmarks/bench_*.py`) are excluded from default
  collection; run with `pytest benchmarks/`. The performance budget is
  "<5 min on M2 chr22-subset"; CI fails on >25% regression.

---

## 11. FAQ (tribal answers)

**Why bundle the chain file?**
The UCSC liftover chains are tiny (223 KB) and rarely change. Bundling
removes a manual install step that's easy to get wrong (UCSC offers
both `hg19ToHg38` and `hg38ToHg19` on the same download page; grabbing
the wrong direction silently breaks Stage 1 yield). SHA verification
catches corruption + tampering.

**Why is the JSON sidecar separate from the JSON report?**
Different consumers. The `<prefix>.pseudohaploid.json` sidecar is
read by **pgen-samplebind at merge time** — it needs to ship alongside
the EIGENSTRAT triplet wherever it goes. The run-report JSON is
**framework telemetry** consumed by ancestry-pipeline-tool's gate
node; it's tied to the run, not the data. Coupling them would force
pgen-samplebind to depend on framework-internal schema.

**Why hg19 as the internal frame?**
AADR is published in hg19 coordinates. The canonical AADR rsID join
is hg19-native; flipping the internal frame to hg38 would require
lifting AADR forward (extra work) and then lifting back at the end
(downstream qpAdm/qpgraph workflows expect hg19 coords).

**Why no CRAM fast path?**
The "fast path" terminology in this codebase means *no-lift*, not
CRAM-related. The no-lift fast path is build-based (hg19 BAM + hg19
AADR), not format-based. A hg38 CRAM still needs Stages 1/2/4 because
the AADR is hg19. CRAM format support is at the pysam layer
(`AlignmentFile` mode `"r"` auto-detects BAM vs CRAM).

**Why the inverted BAQ flag?**
Default behaviour is BAQ disabled (matches pileupCaller's recommended
cmdline; mpileup's `-B` is the disable-BAQ flag). The CLI exposes the
*positive* sense (`--enable-baq` opts in). The double inversion
bridges the two conventions; see `extract_cmd.py:233-242`.

**Why is `samtools -@` not used anymore?**
mpileup is single-threaded — the `-@` flag doesn't exist in any
samtools release we ship against (issue #2). The `--threads` CLI flag
is preserved for back-compat but only affects pileupCaller (which has
its own `--threads`).

**Where is `--quiet`?**
Root group, not subcommand. `pileup-aadr --quiet extract …`, not
`pileup-aadr extract --quiet …`. See §8 Gotchas.

---

## 12. Sibling tools in the chain

`pileup-aadr` is the upstream half of a two-tool chain:

```
BAM ──(pileup-aadr extract)──► EIGENSTRAT triplet + sidecar
                                    │
                                    ▼
                              (pgen-samplebind merge)
                                    │
                                    ▼
                              merged PFILE (.pgen/.pvar/.psam)
                                    │
                                    ▼
                              AdmixTools 2 (qpAdm, qpgraph, etc.)
```

- **pgen-samplebind** (https://github.com/carstenerickson/pgen-samplebind) —
  the downstream merge tool. Reads our `<prefix>.pseudohaploid.json`
  sidecar to skip per-sample heterozygosity recomputation. Its issue #2
  pins the sidecar contract.
- **AdmixTools 2 fork** — install
  `carstenerickson/admixtools@production/v1.0`. Upstream `uqrmaie1/admixtools`
  v2.0.10 doesn't read pgen-samplebind's PFILE output cleanly.

The full chain is what the full-chain f2 invariant test
(`tests/integration/test_lld19_full_chain.py`) defends.

---

## 13. Quick-start checklist for a new contributor

1. Skim §1–4 of this file (architecture + codebase tour + walkthrough).
2. Read [CONTRIBUTING.md](CONTRIBUTING.md) (dev install + how to run tests).
3. Run `pytest` — should collect 313 tests; integration tests auto-skip
   without their toolchain, so it finishes in ~5s either way. If it
   doesn't, fix your env before doing anything else.
4. Browse [open issues](https://github.com/carstenerickson/pileup-aadr/issues)
   for a small starter. (There are no `# TODO:` markers in-tree — the
   project keeps work items in GitHub.)
5. For your first change: trace the call chain in §4 for the area
   you're touching. Write a failing test FIRST, then make it pass.
6. Before committing: run the lint + type-check command from
   [CONTRIBUTING.md § Lint and type-checking](CONTRIBUTING.md#lint-and-type-checking).
   CI will fail otherwise.

Welcome.
