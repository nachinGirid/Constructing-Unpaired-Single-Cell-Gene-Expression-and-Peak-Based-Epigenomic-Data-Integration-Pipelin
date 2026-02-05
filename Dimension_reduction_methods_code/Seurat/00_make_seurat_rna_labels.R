#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: 00_make_seurat_rna_labels.R <rna_counts_rds> <out_label_rds>")
}

inpath  <- args[1]
outpath <- args[2]

rna.data <- readRDS(inpath)
rna <- CreateSeuratObject(counts = rna.data, min.cells = 3)

rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)

rna <- RunPCA(rna, features = VariableFeatures(object = rna))

rna <- FindNeighbors(rna, dims = 1:10)
rna <- FindClusters(rna, resolution = 0.5)

saveRDS(rna@active.ident, outpath)
