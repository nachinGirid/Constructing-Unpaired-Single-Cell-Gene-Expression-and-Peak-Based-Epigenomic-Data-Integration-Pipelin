#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: 00_make_seurat_rna_labels.R <rna_counts_rds> <out_label_rds> <seed> [n_pcs] [resolution]\n",
    "  rna_counts_rds: dgCMatrix genes x cells\n",
    "  out_label_rds: output .rds file for RNA labels (named by cell barcodes)\n",
    "  seed: integer random seed for reproducibility\n",
    "  n_pcs (optional): number of PCs for clustering, default=10\n",
    "  resolution (optional): clustering resolution, default=0.5\n",
    sep = ""
  )
}

rna_path   <- args[1]
out_path   <- args[2]
seed       <- as.integer(args[3])
n_pcs      <- ifelse(length(args) >= 4, as.integer(args[4]), 10L)
resolution <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.5)

set.seed(seed)

rna_counts <- readRDS(rna_path)
if (!inherits(rna_counts, "dgCMatrix") && !is.matrix(rna_counts)) {
  stop("Input RNA counts must be a dgCMatrix (recommended) or a matrix.")
}

rna <- CreateSeuratObject(counts = rna_counts, min.cells = 3)

rna <- NormalizeData(rna, verbose = FALSE)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
rna <- ScaleData(rna, features = rownames(rna), verbose = FALSE)
rna <- RunPCA(rna, features = VariableFeatures(rna), verbose = FALSE)

rna <- FindNeighbors(rna, dims = 1:n_pcs, verbose = FALSE)
rna <- FindClusters(rna, resolution = resolution, verbose = FALSE)

rna_labels <- rna@active.ident
saveRDS(rna_labels, out_path)

cat("Saved RNA labels to: ", out_path, "\n", sep = "")
cat("Clusters: ", paste(levels(rna_labels), collapse = ", "), "\n", sep = "")
