#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Signac)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: 01_preprocess_seurat_objects.R <rna_counts_rds> <atac_counts_rds> <out_dir>\n",
    "  rna_counts_rds: dgCMatrix genes x cells\n",
    "  atac_counts_rds: dgCMatrix peaks x cells\n",
    "  out_dir: output folder\n",
    sep = ""
  )
}

rna_path  <- args[1]
atac_path <- args[2]
out_dir   <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rna_counts  <- readRDS(rna_path)
atac_counts <- readRDS(atac_path)

stopifnot(inherits(rna_counts, "dgCMatrix") || is.matrix(rna_counts))
stopifnot(inherits(atac_counts, "dgCMatrix") || is.matrix(atac_counts))

rna  <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
atac <- CreateSeuratObject(counts = atac_counts, assay = "peaks")

# Optional: remove all-zero features
rna  <- rna[rowSums(rna@assays$RNA@counts) > 0, ]
atac <- atac[rowSums(atac@assays$peaks@counts) > 0, ]

# RNA preprocessing
DefaultAssay(rna) <- "RNA"
rna <- NormalizeData(rna, verbose = FALSE)
rna <- FindVariableFeatures(rna, verbose = FALSE)
rna <- ScaleData(rna, verbose = FALSE)
rna <- RunPCA(rna, verbose = FALSE)

# ATAC preprocessing (TF-IDF + LSI)
DefaultAssay(atac) <- "peaks"
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)  # creates 'lsi'

saveRDS(rna,  file.path(out_dir, "pbmc_preprocessed_seurat_obj_rna.rds"))
saveRDS(atac, file.path(out_dir, "pbmc_preprocessed_seurat_obj_atac.rds"))

cat("Saved:\n",
    " - ", file.path(out_dir, "pbmc_preprocessed_seurat_obj_rna.rds"), "\n",
    " - ", file.path(out_dir, "pbmc_preprocessed_seurat_obj_atac.rds"), "\n", sep = "")
