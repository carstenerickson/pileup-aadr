#!/usr/bin/env bash
# Regenerate the LLD #19 fixtures (AADR ref + pre-extracted triplet + f2 baseline).
#
# Required toolchain on PATH (any Linux box; verified on n2-standard-16):
#   - Python 3.11+ with pileup-aadr installed (any version >= 0.1.2)
#   - pgen-samplebind 0.1.0+
#   - samtools 1.16+, pileupCaller 1.6.0+, picard 3.0+, java 11+
#   - R 4.0+ with admixtools (carstenerickson/admixtools@production/v1.0)
#   - convertf, plink2 (EIGENSOFT + plink2)
#
# Required env vars:
#   PILEUP_AADR_LLD19_BAM   chr22-coverage WGS BAM in hg38 (DRAGEN/bwa/etc.)
#   PILEUP_AADR_LLD19_REF   matching hg38 FASTA (.fai + .dict alongside)
#   PILEUP_AADR_LLD19_AADR  full AADR v66 1240k PACKEDANCESTRYMAP prefix
#   PICARD_JAR              Picard 3.x JAR path
#
# Usage:
#   PILEUP_AADR_LLD19_BAM=/path/to/wgs.bam \
#   PILEUP_AADR_LLD19_REF=/path/to/hg38.fa \
#   PILEUP_AADR_LLD19_AADR=/path/to/v66.1240K.aadr.PUB \
#   PICARD_JAR=/path/to/picard.jar \
#   bash scripts/regen_lld19_baseline.sh
#
# Output: tests/integration/fixtures/lld19/{aadr_ref,pileup_aadr_output,baseline}/*
#
# Wallclock: ~3-5 min on the reference instance (60s for pileup-aadr extract
# + ~30s for AT2 read of full AADR + ~10s for everything else).

set -euo pipefail

: "${PILEUP_AADR_LLD19_BAM:?must be set}"
: "${PILEUP_AADR_LLD19_REF:?must be set}"
: "${PILEUP_AADR_LLD19_AADR:?must be set}"
: "${PICARD_JAR:?must be set}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIX_DIR="$REPO_ROOT/tests/integration/fixtures/lld19"
WORK_DIR="${PILEUP_AADR_LLD19_WORK:-$(mktemp -d -t lld19-regen-XXXX)}"
echo "[regen] work dir: $WORK_DIR"
echo "[regen] fixture dir: $FIX_DIR"

mkdir -p "$FIX_DIR/aadr_ref" "$FIX_DIR/pileup_aadr_output" "$FIX_DIR/baseline"
cd "$WORK_DIR"

# --- Step 1: subset AADR to chr22 + 5 reference populations via AT2 + plink2 ---
echo "[regen] step 1: AADR chr22 5-pop subset"
Rscript -e "
suppressPackageStartupMessages(library(admixtools))
pops <- c('Luxembourg_Loschbour_Mesolithic',
          'Germany_ViesenhaeuserHof_EN',
          'Russia_Malta_UP',
          'Russia_Denisova',
          'Altai_Neanderthal')
packedancestrymap_to_plink(
  inpref='$PILEUP_AADR_LLD19_AADR',
  outpref='aadr_5pops_full',
  pops=pops, verbose=FALSE)
"
plink2 --bfile aadr_5pops_full --chr 22 --make-bed --out aadr_5pops_chr22 \
  --silent --threads 4 >/dev/null

# Convert PLINK to EIGENSTRAT; manually rewrite .ind because convertf can't
# preserve FID→POP mapping cleanly.
cat > convertf.par <<EOF
genotypename:    aadr_5pops_chr22.bed
snpname:         aadr_5pops_chr22.bim
indivname:       aadr_5pops_chr22.fam
genotypeoutname: aadr_chr22_5pops.geno
snpoutname:      aadr_chr22_5pops.snp
indivoutname:    aadr_chr22_5pops.ind.tmp
inputformat:     PACKEDPED
outputformat:    EIGENSTRAT
familynames:     NO
EOF
convertf -p convertf.par > convertf.log 2>&1
awk '{ sex = ($5==1) ? "M" : ($5==2 ? "F" : "U"); print $2, sex, $1 }' \
  aadr_5pops_chr22.fam > aadr_chr22_5pops.ind

cp aadr_chr22_5pops.{geno,snp,ind} "$FIX_DIR/aadr_ref/"
echo "[regen] AADR chr22 ref: $(wc -l < aadr_chr22_5pops.snp) SNPs × $(wc -l < aadr_chr22_5pops.ind) samples"

# --- Step 2: subset BAM to chr22 (skip if already a chr22-only BAM) ---
echo "[regen] step 2: chr22 BAM subset"
samtools view -b -@ 8 "$PILEUP_AADR_LLD19_BAM" chr22 > chr22.bam
samtools index chr22.bam

# --- Step 3: pileup-aadr extract (the canonical 4-stage pipeline) ---
echo "[regen] step 3: pileup-aadr extract"
pileup-aadr extract chr22.bam aadr_chr22_5pops.snp \
  -o extract_out \
  --bam-build hg38 --aadr-build hg19 \
  --ref-fasta "$PILEUP_AADR_LLD19_REF" \
  --picard-mem 8g \
  --report-json extract.report.json \
  --min-coverage 100 --warn-coverage 1000 \
  --liftover-yield-fail-pct 50 \
  --sample-name GFX0442453 --pop GFX0442453_test \
  --quiet > /dev/null

cp extract_out.{geno,snp,ind,pseudohaploid.json} "$FIX_DIR/pileup_aadr_output/"
cp extract.report.json "$FIX_DIR/baseline/"
echo "[regen] pileup-aadr output: $(wc -l < extract_out.geno) variants × 1 sample"

# --- Step 4: pgen-samplebind merge (AADR ref + new sample → merged PFILE) ---
echo "[regen] step 4: pgen-samplebind merge"
pgen-samplebind merge aadr_chr22_5pops --target extract_out -o merged --quiet > /dev/null
echo "[regen] merged: $(wc -l < merged.psam) samples × $(wc -l < merged.pvar) variants (incl. headers)"

# --- Step 5: AT2 extract_f2 → freeze baseline .rds ---
echo "[regen] step 5: AT2 extract_f2"
Rscript -e "
suppressPackageStartupMessages(library(admixtools))
extract_f2('merged', outdir='f2_dir', maxmiss=0, overwrite=TRUE,
           verbose=FALSE, auto_only=FALSE)
f2_blocks <- f2_from_precomp('f2_dir', verbose=FALSE)
saveRDS(f2_blocks, '$FIX_DIR/baseline/f2_baseline.rds')
cat(sprintf('[regen] f2 baseline: %s × %s × %s blocks (%d bytes)\n',
            dim(f2_blocks)[1], dim(f2_blocks)[2], dim(f2_blocks)[3],
            file.size('$FIX_DIR/baseline/f2_baseline.rds')))
"

echo "[regen] done. fixtures regenerated in: $FIX_DIR"
echo "[regen] now run: pytest tests/integration/test_lld19_full_chain.py -v"
