#!/usr/bin/env Rscript

suppressMessages({
  library(Seurat)
  library(Signac)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: cca_transfer.R <dataset> <gas> <dimension>")

ds  <- args[1]
gas <- args[2]
dim <- as.integer(args[3])

basepath <- "processed_data"

# inputs
rna_label <- readRDS(file.path(basepath, ds, "seurat_RNA_label.rds"))
rna  <- readRDS(file.path(basepath, ds, "Seurat_preprocess_rna_obj.rds"))
atac <- readRDS(file.path(basepath, ds, "Seurat_preprocess_atac_obj.rds"))
gas_mat <- readRDS(file.path("GAS/results", sprintf("%s_%s.rds", ds, gas)))

# add GAS
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gas_mat)
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac, verbose = FALSE)
atac <- ScaleData(atac, features = rownames(atac), verbose = FALSE)

# CCA anchors
transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = atac,
  features = VariableFeatures(object = rna),
  reference.assay = "RNA",
  query.assay = "ACTIVITY",
  reduction = "cca",
  dims = 1:dim
)

# label transfer
celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna_label,
  weight.reduction = atac[["lsi"]],
  dims = 2:30
)

# extract CCA embeddings
em <- transfer.anchors@object.list[[1]]@reductions$cca@cell.embeddings
rna_idx  <- grepl("_reference$", rownames(em))
atac_idx <- grepl("_query$", rownames(em))
rna_em  <- em[rna_idx, , drop = FALSE]
atac_em <- em[atac_idx, , drop = FALSE]

# remove suffixes
rownames(rna_em)  <- sub("_reference$", "", rownames(rna_em))
rownames(atac_em) <- sub("_query$", "", rownames(atac_em))

# align by shared barcodes
common <- intersect(rownames(rna_em), rownames(atac_em))
rna_em  <- rna_em[common, , drop = FALSE]
atac_em <- atac_em[common, , drop = FALSE]

# outputs
dir.create("DMR/results", showWarnings = FALSE, recursive = TRUE)
dir.create("CLUSTER/cca_anchor/results", showWarnings = FALSE, recursive = TRUE)

cca <- list(rna = rna_em, atac = atac_em)
saveRDS(cca, sprintf("DMR/results/%s_cca_%s_dim%d.rds", ds, gas, dim))
utils::write.csv(rna_em,  sprintf("DMR/results/%s_rna_%s_dim%d.csv", ds, gas, dim), quote = FALSE, row.names = TRUE)
utils::write.csv(atac_em, sprintf("DMR/results/%s_atac_%s_dim%d.csv", ds, gas, dim), quote = FALSE, row.names = TRUE)

label <- list(
  rna_label = rna_label,
  atac_label = celltype.predictions$predicted.id
)
saveRDS(label, sprintf("CLUSTER/cca_anchor/results/%s_labels_%s_dim%d.rds", ds, gas, dim))
