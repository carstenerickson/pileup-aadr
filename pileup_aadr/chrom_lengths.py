"""Hardcoded hg19 + hg38 chromosome-length tables.

These are byte-stable facts of the assembly versions — they have not changed since the
assemblies were finalized and will not change. Test #30 asserts the tables match UCSC's
published values exactly, catching any future maintainer who accidentally bumps a number.
"""
from typing import Final

HG19_CHROM_LENGTHS: Final[dict[str, int]] = {
    "chr1": 249_250_621,
    "chr2": 243_199_373,
    "chr3": 198_022_430,
    "chr4": 191_154_276,
    "chr5": 180_915_260,
    "chr6": 171_115_067,
    "chr7": 159_138_663,
    "chr8": 146_364_022,
    "chr9": 141_213_431,
    "chr10": 135_534_747,
    "chr11": 135_006_516,
    "chr12": 133_851_895,
    "chr13": 115_169_878,
    "chr14": 107_349_540,
    "chr15": 102_531_392,
    "chr16": 90_354_753,
    "chr17": 81_195_210,
    "chr18": 78_077_248,
    "chr19": 59_128_983,
    "chr20": 63_025_520,
    "chr21": 48_129_895,
    "chr22": 51_304_566,
    "chrX": 155_270_560,
    "chrY": 59_373_566,
    "chrM": 16_571,
}

HG38_CHROM_LENGTHS: Final[dict[str, int]] = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_441,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
    "chrX": 156_040_895,
    "chrY": 57_227_415,
    "chrM": 16_569,
}

# Iteration order for ##contig headers (matches Picard's preference: numeric first,
# then X/Y/M; alphabetical within numerics is wrong because chr11 sorts before chr2).
CHROM_ORDER: Final[list[str]] = (
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
)

__all__ = ["CHROM_ORDER", "HG19_CHROM_LENGTHS", "HG38_CHROM_LENGTHS"]
