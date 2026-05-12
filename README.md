# pileup-aadr

> One-shot user-BAM → AADR-site pseudohaploid genotypes for ancient-DNA pipelines.

A focused tool that takes a user's WGS BAM (hg19 or hg38) plus an AADR `.snp` file and emits coverage-matched pseudohaploid genotypes in EIGENSTRAT format, ready for sample-binding via [pgen-samplebind](https://github.com/carstenerickson/pgen-samplebind). Replaces the 5-step Picard-liftover + sites-VCF-roundtrip + pileupCaller + rsID-rejoin dance that every ancient-DNA personal-WGS pipeline currently reimplements.

**Status: alpha, in active development.** Day 1 of 15 of the v0.1.0 implementation. See [Project plan](#project-plan) below.

## Why it exists

The ancient-DNA community has converged on AADR's 1240k SNP capture panel + EIGENSOFT `convertf` + Reich Lab's `pileupCaller` + `mergeit` as the canonical pipeline for getting personal/cohort WGS data into the qpAdm/qpgraph workflow. Each user reimplements the same 5-step process:

1. Lift AADR `.snp` (hg19) → sites VCF in user's build (hg38) via Picard
2. Convert sites VCF → pileupCaller `.snp` format
3. Run `samtools mpileup` against the user's BAM at lifted sites
4. Run `pileupCaller --randomDiploid` to call pseudohaploid genotypes
5. Rejoin to hg19 coordinates via rsID for AdmixTools 2 compatibility

The procedure is documented in scattered methodology notes; execution is error-prone (~99% liftover yield is required; lower numbers indicate chain or reference-FASTA mismatch); coverage reporting is ad-hoc. Failure modes ("only 39% of HO sites lifted cleanly" via the gVCF-lift dead-end) are tribal knowledge.

`pileup-aadr` is a single-purpose tool wrap of these 5 steps with proper instrumentation + validation, removing the most failure-prone manual step in any ancient-DNA pipeline.

## Subcommands

```
pileup-aadr extract  BAM AADR_SNP -o PREFIX [...]   the canonical 4-stage pipeline
pileup-aadr validate BAM AADR_SNP                   pre-flight check (no mpileup/pileupCaller)
pileup-aadr coverage BAM [--regions BED]            per-chrom coverage report (mosdepth)
pileup-aadr inspect  AADR_SNP                       structured summary of an AADR .snp panel
```

## Install (alpha)

```bash
pip install -e ".[dev]"   # editable install with dev deps
```

External binaries needed by `extract`:
- `samtools` ≥ 1.16
- `pileupCaller` ≥ 1.6.0
- `picard` ≥ 3.0 (only when AADR build ≠ BAM build; not needed for hg19-native fast path)
- `java` ≥ 11 (transitive Picard requirement)

For `coverage` subcommand:
- `mosdepth` ≥ 0.3.6

Easiest install (most CI matrices): `conda install -c bioconda samtools pileupcaller picard mosdepth`.

## Design

This repo implements a frozen high-level + low-level design (kept in private notes).
Key invariants the implementation honors:

- **4-stage extract pipeline** mirroring the methodology's Step 1.4–1.7 (lift → transform → call → rejoin)
- **Picard 3.3.0+ LiftoverVcf** with `RECOVER_SWAPPED_REF_ALT=true` + `SwappedAlleles` INFO flag for Stage 4 dosage inversion
- **Bundled UCSC chain** (~223 KB; SHA-verified at startup)
- **`<prefix>.pseudohaploid.json` sidecar** for authoritative provenance (consumed by `pgen-samplebind`'s sidecar reader; see [pgen-samplebind#2](https://github.com/carstenerickson/pgen-samplebind/issues/2))
- **Stable exit codes** (0/1/2/3/4) for workflow-manager integration

CHANGELOG entries summarize the per-day implementation progress + reference the
relevant design constraints inline.

## Project plan (v0.1.0)

3 weeks; 15 working days.

**Week 1**: skeleton + format detection + Stage 1 (Picard lift) + Stage 2 (transform) + Stage 3 (mpileup|pileupCaller pipe) + Stage 4 (rsID rejoin) + correctness invariants.

**Week 2**: hg19-native fast path + ancestry-pipeline-tool integration (stdout/stderr discipline, exit-code map, JSON schema) + pgen-samplebind composition test + robustness regression + concurrency.

**Week 3**: packaging + self-dogfood + external feedback iteration + v0.1.0 tag.

CHANGELOG tracks day-by-day progress.

## Status

- Design: frozen, all readiness items closed (sidecar integration confirmed via [pgen-samplebind#2](https://github.com/carstenerickson/pgen-samplebind/issues/2))
- Implementation: Day 1 in progress

## License

MIT — see [LICENSE](LICENSE).
