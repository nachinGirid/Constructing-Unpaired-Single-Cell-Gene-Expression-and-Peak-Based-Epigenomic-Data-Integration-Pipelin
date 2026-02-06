#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(bindSC)
  library(Matrix)
})

### For example pbmc10k data
ds  <- "pbmc"
gas <- "signac"

# input files (user should edit paths)
rna  <- readRDS("seurat_rna_obj.rds")
atac <- readRDS("seurat_atac_obj.rds")
gasm <- readRDS(sprintf("%s_%s.rds", ds, gas))   # genes x cells (dgCMatrix)

# -------- align GAS to ATAC cells --------
common_cells <- intersect(colnames(gasm), colnames(atac))
gasm <- gasm[, common_cells, drop = FALSE]
atac <- subset(atac, cells = common_cells)
rna  <- subset(rna,  cells = common_cells)

# -------- RNA: normalize + HVGs --------
DefaultAssay(rna) <- "RNA"
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = 5000)

# -------- add ACTIVITY assay and normalize --------
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gasm)
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- FindVariableFeatures(atac, nfeatures = 5000)

gene.use <- intersect(
  VariableFeatures(rna),
  VariableFeatures(atac)
)

# -------- build X, Z0 --------
X  <- GetAssayData(rna,  assay = "RNA",      slot = "data")[gene.use, ]
Z0 <- GetAssayData(atac, assay = "ACTIVITY", slot = "data")[gene.use, ]

# -------- peaks LSI for Y --------
DefaultAssay(atac) <- "ATAC"
atac <- RunTFIDF(atac)
acc  <- Matrix::rowSums(GetAssayData(atac, assay = "ATAC", slot = "counts"))
keep <- which(acc > 50)
atac <- RunSVD(atac, n = 50, features = rownames(atac)[keep])
y <- atac@reductions$lsi@cell.embeddings

# -------- joint pre-reduction --------
dr <- dimReduce(dt1 = X, dt2 = Z0, K = 30)
x  <- dr$dt1
z0 <- dr$dt2

# -------- pre-cluster (BiCCA requires) --------
# ATAC clustering
atac <- FindNeighbors(atac, reduction = "lsi", dims = 1:20)
atac <- FindClusters(atac, resolution = 0.5)
y.clst <- atac$seurat_clusters

# RNA clustering
rna1 <- CreateSeuratObject(counts = X)
rna1 <- NormalizeData(rna1)
rna1 <- FindVariableFeatures(rna1, nfeatures = 5000)
rna1 <- ScaleData(rna1, features = rownames(rna1))
rna1 <- RunPCA(rna1, npcs = 50)
rna1 <- FindNeighbors(rna1, dims = 1:20)
rna1 <- FindClusters(rna1, resolution = 0.5)
x.clst <- rna1$seurat_clusters

# -------- BiCCA --------
res <- BiCCA(
  X = t(x),
  Y = t(y),
  Z0 = t(z0),
  X.clst = x.clst,
  Y.clst = y.clst,
  alpha = 0.5,
  lambda = 0.5,
  K = 15,
  num.iteration = 50,
  tolerance = 0.01,
  save = TRUE,
  parameter.optimize = FALSE,
  block.size = 0
)

U <- res$u
R <- res$r

if (is.null(colnames(U))) colnames(U) <- paste0("Dim", seq_len(ncol(U)))
if (is.null(colnames(R))) colnames(R) <- paste0("Dim", seq_len(ncol(R)))

# outputs
utils::write.csv(
  U,
  sprintf("%s_rna_%s_RNA_embedding.csv", ds, gas),
  quote = FALSE,
  row.names = TRUE
)

utils::write.csv(
  R,
  sprintf("%s_atac_%s_ATAC_embedding.csv", ds, gas),
  quote = FALSE,
  row.names = TRUE
)
