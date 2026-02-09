#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

### For example pbmc data
method <- "seurat"
ds     <- "pbmc"

# ----------------------------
# User-editable paths
# ----------------------------
basepath  <- "PATH/TO/processed_embeddings"
data_root <- "PATH/TO/processed_data"

# FigR helper functions (user should edit paths)
source("cellPairing.R")
source("utils.R")
source("FigR.R")
source("Function_anchor_transfer_label.r")

# one-shot auto knob-----------------------------------------------------------------------------
choose_figr_knobs <- function(n_atac, n_rna) {
  N <- (n_atac + n_rna)/2

  # search_range controls KNN window size (search_range × total cells)
  search_range <- if (N <= 5000) {
    0.2
  } else if (N <= 30000) {
    0.05
  } else {
    0.005
  }

  # tune as need
  max_multimatch <- 10L #30, 50

  # minimum subgraph size scales with dataset size
  min_subgraph_size <- if (N <= 5000) {
    50L
  } else if (N <= 30000) {
    100L
  } else {
    200L
  }

  list(
    search_range      = search_range,
    max_multimatch    = max_multimatch,
    min_subgraph_size = min_subgraph_size
  )
}

# ----------------------------
# Input files
# ----------------------------
embed_dir <- file.path(basepath, method, ds)

rna_p  <- file.path(embed_dir, "rna_em.csv")
atac_p <- file.path(embed_dir, "atac_em.csv")

# output
out_p <- file.path(embed_dir, "atac_labels_figR.txt")

# ----------------------------
# Helpers
# ----------------------------
read_mat <- function(p) {
  df <- fread(p, data.table = FALSE, check.names = FALSE)
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  as.matrix(df)
}

read_rna_label_with_barcodes <- function(ds) {
  lab_txt <- file.path(data_root, ds, "seurat_RNA_label.txt")
  lab_csv <- file.path(data_root, ds, "seurat_RNA_label.csv")
  bc_p    <- file.path(data_root, ds, "barcodes.tsv")

  lab_p <- if (file.exists(lab_txt)) lab_txt else lab_csv
  lab <- fread(lab_p, data.table = FALSE)
  if (ncol(lab) != 1) stop("RNA label file must have one column")

  bc <- fread(bc_p, header = FALSE, data.table = FALSE)[, 1]
  rownames(lab) <- bc
  lab[[1]]
}

# ----------------------------
# Main
# ----------------------------
rna_em  <- read_mat(rna_p)
atac_em <- read_mat(atac_p)
rna_lab <- read_rna_label_with_barcodes(ds)

common <- Reduce(intersect, list(
  rownames(rna_em),
  rownames(atac_em),
  names(rna_lab)
))

rna_em  <- rna_em[common, , drop = FALSE]
atac_em <- atac_em[common, , drop = FALSE]
rna_lab <- rna_lab[common]

# suffix modality
rownames(rna_em)  <- paste0(rownames(rna_em), "_rna")
rownames(atac_em) <- paste0(rownames(atac_em), "_atac")

# FigR pairing
pairing <- pairCells(
  ATAC              = atac_em,
  RNA               = rna_em,
  search_range      = search_range,
  max_multimatch    = max_multimatch,
  keepUnique        = TRUE,
  min_subgraph_size = min_subgraph_size
)

# label transfer
atac_label <- transfer_rna_labels_to_atac(
  pairing          = pairing,
  rna_lab          = rna_lab,
  atac_pcs         = atac_em,
  neighbors_number = neighbors_number
)

# save
out_df <- data.frame(
  label = atac_label,
  row.names = names(atac_label),
  check.names = FALSE
)

fwrite(
  as.data.table(setDT(out_df, keep.rownames = TRUE)),
  file = out_p,
  sep = "\t",
  col.names = TRUE
)

cat("Done: FigR label transfer\n")
