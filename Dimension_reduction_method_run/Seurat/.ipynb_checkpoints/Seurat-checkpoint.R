suppressMessages({
  library(Seurat)
  library(Signac)
  library(stringr)
})
### For example pbmc10k data
ds  <- 'pbmc'
gas <- 'signac'
dim <- 15

basepath <- "PATH/TO/PROCESSED_DATA"   # user should edit

rna_label <- readRDS(file.path(basepath, ds, "seurat_RNA_label.rds"))
rna       <- readRDS(file.path(basepath, ds, "Seurat_preprocess_rna_obj.rds"))
atac      <- readRDS(file.path(basepath, ds, "Seurat_preprocess_atac_obj.rds"))

# GAS matrix (method-specific)
gas_mat   <- readRDS(sprintf("%s_%s.rds", ds, gas))

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

# label transfer (uses precomputed LSI in atac if present)
celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna_label,
  weight.reduction = atac[["lsi"]],
  dims = 2:30
)

# extract CCA embeddings
em <- transfer.anchors@object.list[[1]]@reductions$cca@cell.embeddings

rna_idx  <- grepl("_reference$", rownames(em))
atac_idx <- grepl("_query$",     rownames(em))

rna_em  <- em[rna_idx,  , drop = FALSE]
atac_em <- em[atac_idx, , drop = FALSE]

# remove suffixes
rownames(rna_em)  <- sub("_reference$", "", rownames(rna_em))
rownames(atac_em) <- sub("_query$",     "", rownames(atac_em))

# align by shared barcodes
common <- intersect(rownames(rna_em), rownames(atac_em))
rna_em  <- rna_em[common, , drop = FALSE]
atac_em <- atac_em[common, , drop = FALSE]

utils::write.csv(rna_em,sprintf("%s_rna_%s_dim%d.csv", ds, gas, dim),quote = FALSE, row.names = TRUE)
utils::write.csv(atac_em,sprintf("%s_atac_%s_dim%d.csv", ds, gas, dim),quote = FALSE, row.names = TRUE)


## we also save the seurat -anchor based label transferring result from this script
label <- list(
  rna_label = rna_label,
  atac_label = celltype.predictions$predicted.id
)
saveRDS(label,sprintf("%s_labels_%s_dim%d.rds", ds, gas, dim))
