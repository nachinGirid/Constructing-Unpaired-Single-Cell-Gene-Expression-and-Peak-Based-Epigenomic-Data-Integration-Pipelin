#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 01_preprocess_seurat_objects.R <rna_counts_rds> <atac_counts_rds> <out_dir>")
}

rna_inpath  <- args[1]
atac_inpath <- args[2]
out_dir     <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rna.data  <- readRDS(rna_inpath)
atac.data <- readRDS(atac_inpath)

rna  <- CreateSeuratObject(counts = rna.data)
atac <- CreateSeuratObject(counts = atac.data)

saveRDS(rna,  file.path(out_dir, "seurat_rna_obj.rds"))
saveRDS(atac, file.path(out_dir, "seurat_atac_obj.rds"))

atac_fil <- atac[rowSums(atac) != 0, ]
rna_fil  <- rna[rowSums(rna) != 0, ]

saveRDS(rna_fil,  file.path(out_dir, "Seurat_obj_rna_fil_rowsum0.rds"))
saveRDS(atac_fil, file.path(out_dir, "Seurat_obj_atac_fil_rowsum0.rds"))

# RNA preprocess
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna)
rna <- ScaleData(rna)
rna <- RunPCA(rna)

# ATAC preprocess
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)

saveRDS(rna,  file.path(out_dir, "Seurat_preprocess_rna_obj.rds"))
saveRDS(atac, file.path(out_dir, "Seurat_preprocess_atac_obj.rds"))
