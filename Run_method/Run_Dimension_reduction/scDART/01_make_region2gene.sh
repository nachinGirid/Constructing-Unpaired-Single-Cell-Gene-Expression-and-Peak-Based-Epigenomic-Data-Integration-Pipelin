#!/usr/bin/env bash

### For example pbmc10k data

# input files (user should edit paths)
GENE_BED="hg38_genebody_upstream2000.bed"
PEAK_BED="pbmc_peak.bed"

# output
OUT="pbmc_region2gene.txt"

bedtools intersect \
  -a ${GENE_BED} \
  -b ${PEAK_BED} \
  -wa -wb | \
  cut -f 5-8 > ${OUT}
