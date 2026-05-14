# pileup-aadr

> One-shot user-BAM → AADR-site pseudohaploid genotypes for ancient-DNA pipelines.

A focused tool that takes a user's WGS BAM (hg19 or hg38) plus an AADR `.snp` file and emits coverage-matched pseudohaploid genotypes in EIGENSTRAT format, ready for sample-binding via [pgen-samplebind](https://github.com/carstenerickson/pgen-samplebind). Replaces the 5-step Picard-liftover + sites-VCF-roundtrip + pileupCaller + rsID-rejoin dance that every ancient-DNA personal-WGS pipeline currently reimplements.

**Status: alpha (v0.1 surface complete).** All four subcommands functional; 208 unit + integration tests passing on a 6-cell CI matrix (Python 3.11/3.12/3.13 × Linux/macOS).

## Why it exists

The ancient-DNA community has converged on AADR's 1240k SNP capture panel + EIGENSOFT `convertf` + Reich Lab's `pileupCaller` + `mergeit` as the canonical pipeline for getting personal/cohort WGS data into the qpAdm/qpgraph workflow. Each user reimplements the same 5-step process:

1. Lift AADR `.snp` (hg19) → sites VCF in user's build (hg38) via Picard
2. Convert sites VCF → pileupCaller `.snp` format
3. Run `samtools mpileup` against the user's BAM at lifted sites
4. Run `pileupCaller --randomDiploid` to call pseudohaploid genotypes
5. Rejoin to hg19 coordinates via rsID for AdmixTools 2 compatibility

The procedure is documented in scattered methodology notes; execution is error-prone (~99% liftover yield is required; lower numbers indicate chain or reference-FASTA mismatch); coverage reporting is ad-hoc. Failure modes ("only 39% of HO sites lifted cleanly" via the gVCF-lift dead-end) are tribal knowledge.

`pileup-aadr` wraps these 5 steps with proper instrumentation + validation, removing the most failure-prone manual step in any ancient-DNA pipeline.

## Install

```bash
pip install -e ".[dev]"
```

External binaries needed by `extract`:

| Tool          | Min version | Required for                                    |
|---------------|-------------|--------------------------------------------------|
| `samtools`    | ≥ 1.16      | always (Stage 3 mpileup)                         |
| `pileupCaller`| ≥ 1.6.0     | always (Stage 3 randomDiploid caller)           |
| `picard`      | ≥ 3.0       | only when AADR build ≠ BAM build (Stage 1 lift)  |
| `java`        | ≥ 11        | transitive Picard requirement                    |
| `mosdepth`    | ≥ 0.3.6     | only `coverage` subcommand                       |

Easiest install: `conda install -c bioconda samtools pileupcaller picard mosdepth`. The UCSC `hg19ToHg38.over.chain.gz` (~223 KB) is bundled with the package and SHA-verified at startup — no separate download needed.

## Quickstart

```bash
# Sanity-check inputs before a 30-min run
pileup-aadr validate /data/sample.bam /data/aadr_v66.snp

# Inspect an AADR panel (panel-size guess, allele distribution, palindrome %)
pileup-aadr inspect /data/aadr_v66.snp

# The main pipeline: BAM + AADR .snp → EIGENSTRAT triplet
pileup-aadr extract /data/sample.bam /data/aadr_v66.snp \
    -o /data/sample_pseudohaploid \
    --report-json /data/sample.report.json

# When the coverage gate fails, triage with mosdepth
pileup-aadr coverage /data/sample.bam --json
```

The `extract` subcommand writes 4 output files at the prefix:

| File                         | Format                          | Purpose                                                                    |
|------------------------------|---------------------------------|----------------------------------------------------------------------------|
| `<prefix>.geno`              | EIGENSTRAT 1-char-per-line      | per-variant pseudohaploid dosage (0/1/2/9)                                  |
| `<prefix>.snp`               | EIGENSTRAT 6-col TSV            | AADR rsID + hg19 chrom + Morgans + pos + REF + ALT                         |
| `<prefix>.ind`               | EIGENSTRAT 3-col TSV            | sample IID + SEX + POP                                                      |
| `<prefix>.pseudohaploid.json`| schema-versioned JSON           | provenance sidecar (consumed by `pgen-samplebind`'s sidecar reader)        |

## Subcommands

```
pileup-aadr extract  BAM AADR_SNP -o PREFIX [...]   the canonical 4-stage pipeline
pileup-aadr validate BAM AADR_SNP                   pre-flight check (no mpileup)
pileup-aadr coverage BAM [--regions BED]            per-chrom coverage report (mosdepth)
pileup-aadr inspect  AADR_SNP                       structured AADR .snp summary
```

Run `pileup-aadr <cmd> --help` for the full option list. Notable `extract` options:

- **Build detection**: `--bam-build` / `--aadr-build` override the auto-detection
- **Output control**: `--sample-name`, `--pop`, `--sex`, `--overwrite`
- **Liftover tuning**: `--chain`, `--ref-fasta`, `--picard-mem` (default 3g), `--strict-chain-sha`
- **Filtering**: `--keep-palindromes`, `--keep-alt-contigs` (rarely useful; defaults are right)
- **Pileup**: `--threads`, `--min-mapq`, `--min-baseq`, `--enable-baq`, `--seed`
- **Gates**: `--liftover-yield-fail-pct` (default 70), `--min-coverage` (default 500k)
- **Reporting**: `--report-json`, `--report-tsv`
- **Tempdir**: `--tempdir`, `--keep-tempdir`, `--clean-tempdir-on-crash`

## Panels: 1240k, HO, chrY/chrM, custom

`pileup-aadr` is panel-agnostic — any EIGENSTRAT `.snp` works. The
`--aadr-snp` argument is just a positional path. `pileup-aadr inspect`
reports a `chrom_set` field that classifies the panel and (for `extract`)
drives whether the autosomal coverage gate applies:

| `chrom_set`      | Typical panel                              | `extract` coverage gate |
|------------------|--------------------------------------------|--------------------------|
| `autosomes+sex`  | AADR v66 1240k, AADR v66 HO                | applies (gate on autosomal calls) |
| `autosomes_only` | 1240k autosomes-only subset                | applies |
| `chrY_only`      | chrY haplogroup workflows                  | **skipped** with INFO log; per-chrom counts in JSON report for downstream sanity |
| `chrM_only`      | mtDNA-only ancient workflows               | **skipped**, same |
| `sex_only`       | chrX + chrY only                           | **skipped**, same |
| `custom`         | one-chrom slices, novel panels             | applies (autosomal threshold; user overrides via `--min-coverage` if needed) |

Tested against AADR v66 1240k (~1.2M sites). HO (~600K sites) shares the
format and works identically; the user is responsible for choosing a
panel-appropriate `--min-coverage` (the 500k default is calibrated for
1240k). For custom panels, run `pileup-aadr inspect <panel.snp>` first
to confirm the `chrom_set` and panel size before pinning thresholds.

For a chrY-only haplogroup workflow:

```bash
# AADR slice with only chr24 (chrY in EIGENSTRAT encoding) rows.
pileup-aadr inspect chrY_panel.snp --json | grep chrom_set
# → "chrom_set": "chrY_only"

pileup-aadr extract sample.bam chrY_panel.snp -o out
# → INFO: Skipping autosomal coverage gate: AADR panel is 'chrY_only' (no
#   autosomes; --min-coverage threshold doesn't apply). Per-chrom call
#   counts available in coverage.per_chrom_call_count for downstream sanity.
```

## Configuration

Two environment variables let you point at a shared site-wide install of the chain + reference FASTAs:

| Env var                    | Effect                                                                |
|----------------------------|------------------------------------------------------------------------|
| `PILEUP_AADR_CHAIN_DIR`    | Directory containing `hg19ToHg38.over.chain.gz` (overrides bundled)   |
| `PILEUP_AADR_REF_DIR`      | Directory containing `<build>.fa` (e.g., `hg38.fa` + `hg38.fa.fai`)   |
| `PILEUP_AADR_JSON_LOGS=1`  | Switch stderr logging from human-readable to JSON Lines                |

Resolution order for `--chain` / `--ref-fasta`: explicit CLI flag → env var → bundled (chain) or BAM @PG (FASTA). Pre-flight verifies the FASTA's chr1 length matches the BAM's build before Picard burns 30+ seconds rejecting most sites with `MismatchedRefAllele`.

## Exit codes

Stable across versions for workflow-manager integration:

| Code | Meaning                                                                 |
|------|--------------------------------------------------------------------------|
| 0    | Success (possibly with `WARN`-level gate flags)                          |
| 1    | Soft-validation failure: liftover yield gate, coverage gate              |
| 2    | I/O failure: chain/FASTA/BAM not found, subprocess crashed, lock held    |
| 3    | Invariant violation: build mismatch, AADR malformed, defensive sanity   |
| 4    | Usage error: bad CLI args, missing/wrong-version external binary         |

## Design

The 4-stage pipeline:

```
AADR .snp (hg19)              user BAM (hg19 or hg38)
       │                            │
       │ Stage 1: Picard LiftoverVcf RECOVER_SWAPPED_REF_ALT
       │ (skipped when BAM build == AADR build)
       ▼
   lifted VCF + SwappedAlleles INFO flags
       │
       │ Stage 2: alt-contig filter + pileupCaller .snp + BED
       ▼
   .snp + BED in target build
       │
       │ Stage 3: samtools mpileup | pileupCaller --randomDiploid
       │ (~25-40 min on a 33× WGS at 1240k)
       ▼
   pileupCaller EIGENSTRAT triplet (target build)
       │
       │ Stage 4: rejoin to AADR's hg19 frame by rsID
       │          + invert dosage at SwappedAlleles flag
       ▼
   final EIGENSTRAT triplet + PSEUDOHAPLOID sidecar
```

Key invariants the implementation honors:

- **Picard 3.3.0+ LiftoverVcf** with `RECOVER_SWAPPED_REF_ALT=true` + `SwappedAlleles` INFO flag for Stage 4 dosage inversion (the safe, modern path; the gVCF-lift route fails at chain-boundary straddling)
- **Bundled UCSC chain** (~223 KB; SHA-verified at startup; reinstall the wheel if the SHA mismatches)
- **Single-BAM `--randomDiploid` is pseudohaploid by construction** — recorded in the `<prefix>.pseudohaploid.json` sidecar for downstream tooling
- **No-lift fast path** when BAM build matches AADR build (saves Picard + transform + rejoin; pileupCaller's output IS the AADR-frame triplet)
- **Stable exit codes** (0/1/2/3/4) for workflow-manager integration
- **Streaming Stage-4 writes** — 1.2M-site EIGENSTRAT triplet emitted at constant memory

## Integration with ancestry-pipeline-tool

`pileup-aadr extract --report-json <path>` produces a schema-versioned JSON consumed by ancestry-pipeline-tool's gate node. The schema (v1) places per-stage `ExtractCounters` fields at the top level alongside `tool`, `input`, `output`, `gates`, and `config` blocks; the no-lift fast path serializes Stages 1/2/4 as `null` so consumers can branch cleanly.

The `<prefix>.pseudohaploid.json` sidecar is read by `pgen-samplebind` via [pgen-samplebind#2](https://github.com/carstenerickson/pgen-samplebind/issues/2) — `pseudohaploid=1` + `het_count` + `het_rate` flow through to the sample-bind step without re-derivation.

## Status

| Surface           | State                                                       |
|-------------------|-------------------------------------------------------------|
| Design            | Frozen (HLD + LLD reconciled; all readiness items closed)   |
| `extract`         | Functional end-to-end (lift + no-lift fast path)            |
| `validate`        | Functional (10-check pre-flight)                            |
| `coverage`        | Functional (mosdepth wrapper)                               |
| `inspect`         | Functional (pure-Python AADR `.snp` summary)                |
| Tests             | 208 passing across 6-cell matrix; ruff-clean                |
| First tag         | Pending real-binary smoke test against captured baselines   |

CHANGELOG tracks the day-by-day implementation progress with design-constraint references.

## Contributing

Dev install, integration-test setup, lint/typing, and the release process
all live in [CONTRIBUTING.md](CONTRIBUTING.md). For a tour of the codebase
itself — architecture, key contracts, common-task recipes — see
[DEVELOPMENT.md](DEVELOPMENT.md). Issues and PRs welcome.

## License

MIT — see [LICENSE](LICENSE).
