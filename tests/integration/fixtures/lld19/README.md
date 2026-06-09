# LLD test #19 fixtures

Fixtures for `tests/integration/test_lld19_full_chain.py` — the
end-to-end correctness invariant from LLD §"pgen-samplebind composition":

> **`test_pileup_aadr_to_samplebind`**: output of `pileup-aadr extract` feeds
> directly into `pgen-samplebind merge --target` without conversion; AT2
> `extract_f2` on the merged panel matches a baseline within `max_dev < 1e-9`.

## Files

| Path | Bytes | Source |
|---|---|---|
| `aadr_ref/aadr_chr22_5pops.{geno,snp,ind}` | 1.2 MB | AADR v66 1240k subset to chr22 + 5 ancient populations (Loschbour, Stuttgart, MA1, Denisova, Altai_Neanderthal). Built via `admixtools::packedancestrymap_to_plink` → `plink2 --chr 22` → `convertf` PACKEDPED → EIGENSTRAT. |
| `pileup_aadr_output/extract_out.{geno,snp,ind,pseudohaploid.json}` | 587 KB | Output of `pileup-aadr extract` (default `--calling-mode randomHaploid`, v0.5.0) against a chr22 subset of an hg38 WGS BAM (DRAGEN-aligned, 33× coverage). 16,188 variants × 1 sample (16,171 non-missing, 0 het — pseudo-haploid by construction). Pre-extracted so the layer-A test doesn't need a BAM. |
| `baseline/f2_baseline.rds` | 2.0 KB | The frozen f2 set: 6×6×13 array (5 ref pops + the new sample × 13 jackknife blocks), produced by AT2 `extract_f2` on the merged panel. The test asserts re-running the chain produces a matching f2 within `max_dev < 1e-9`. |
| `baseline/extract.report.json` | 3.4 KB | Provenance: pileup-aadr's `--report-json` from the run that produced `extract_out.*`. Records counters, gates, tool versions, exit timing. |

## How to regenerate

When pileup-aadr's pipeline changes in a way that legitimately moves the f2 numbers
(e.g., a new Stage 4 invariant), regenerate the baseline via
`scripts/regen_lld19_baseline.sh`. That script encodes the procedure end-to-end on
a host with all four toolchains: AT2 fork (`carstenerickson/admixtools@production/v1.0`),
`pgen-samplebind`, real binaries (samtools/picard/pileupCaller), and a chr22 BAM.

The reference instance for regeneration is `ancestrytracke-f` (n2-standard-16).
Local M-class hardware works too — install AT2 + pgen-samplebind + the bio tools,
point the script at a chr22-subset BAM, and the same fixtures + baseline regenerate
in ~2 minutes.
