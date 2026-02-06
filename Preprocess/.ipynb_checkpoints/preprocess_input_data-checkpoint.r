library(dplyr)
library(Seurat)
library(Signac)
library(SingleCellExperiment)
library(zellkonverter)

# ---- load raw file ----
counts<-Read10X_h5("filtered_feature_bc_matrix.h5")

####### or ####
#counts<-ReadMtx(
#  mtx="matrix.mtx",
#  cells="barcodes.tsv",
#  features="features.tsv")
##################

rna=counts$'Gene Expression'
atac=counts$'Peaks'
# ---- save raw dgCMatrix ----
saveRDS(rna, "rna_raw_dgCMatrix.rds"))
saveRDS(atac, "atac_raw_dgCMatrix.rds"))

# ---- load raw dgCMatrix ----
rna_mat  <- readRDS("rna_raw_dgCMatrix.rds")  
atac_mat <- readRDS("atac_raw_dgCMatrix.rds") 

# ---- Seurat objects ----
rna  <- CreateSeuratObject(counts = rna_mat, min.cells = 3)
atac <- CreateSeuratObject(counts = atac_mat, assay = "ATAC", min.cells = 3)

# ---- save raw Seurat ----
saveRDS(rna, "seurat_rna_obj.rds")
saveRDS(atac, "seurat_atac_obj.rds")


# ---- SCE + H5AD ----
# ---- filter rows with zero counts ----
rna_keep  <- Matrix::rowSums(GetAssayData(rna,  slot = "counts"))  > 0
atac_keep <- Matrix::rowSums(GetAssayData(atac, slot = "counts")) > 0
rna_fil   <- subset(rna,  features = rownames(rna)[rna_keep])
atac_fil  <- subset(atac, features = rownames(atac)[atac_keep])
RNA_sce  <- as.SingleCellExperiment(rna_fil)
ATAC_sce <- as.SingleCellExperiment(atac_fil)
zellkonverter::writeH5AD(RNA_sce, "rna_fil_rowsum0.h5ad")
zellkonverter::writeH5AD(ATAC_sce, "atac_fil_rowsum0.h5ad")

# ---- independent preprocessing ----
# RNA
rna <- NormalizeData(rna, verbose = FALSE)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
rna <- ScaleData(rna, features = rownames(rna), verbose = FALSE)
rna <- RunPCA(rna, features = VariableFeatures(rna), verbose = FALSE)

# ATAC
DefaultAssay(atac) <- "ATAC"
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)

saveRDS(rna, "Seurat_preprocess_rna_obj.rds")
saveRDS(atac,  "Seurat_preprocess_atac_obj.rds")



### some method needs Gene activity score (GAS) as h5ad, we show exampl of a pbmc signac GAS convert it to h5ad:
suppressMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(Matrix)
})

ds  <- "pbmc"
gas <- "signac"

gas_mat <- readRDS("pbmc_signac.rds")  # dgCMatrix genes x cells

atac <- CreateSeuratObject(counts = gas_mat, assay = "ATAC")

# filter zero-row genes (required by many tools)
cnt <- GetAssayData(atac, assay = "ATAC", slot = "counts")
keep_feat <- Matrix::rowSums(cnt) != 0
if (any(!keep_feat)) {
  atac <- subset(atac, features = rownames(atac)[keep_feat])
}

ATAC_sce <- as.SingleCellExperiment(atac)

# write .h5ad
zellkonverter::writeH5AD(ATAC_sce,sprintf("%s_%s.h5ad", ds, gas))
