PCA_ASW <- function(GAS, true_label, n_pcs = 15, scale_factor = 1e4) {
  # GAS: genes x cells
  stopifnot(length(true_label) == ncol(GAS))
  suppressPackageStartupMessages(library(cluster))
  # Keep factor mapping consistent
  labs_factor <- factor(true_label)
  labs  <- as.integer(labs_factor)       # 1..K
  lvls  <- levels(labs_factor)           # cell-type names (in same order)
  
  # 1) Filter genes with zero variance
  sds <- apply(GAS, 1, stats::sd)
  keep <- !is.na(sds) & sds > 0
  gas_fil <- GAS[keep, , drop = FALSE]
  
  # 2) seurat standard
  library(Seurat)
  obj<- CreateSeuratObject(counts = gas_fil)
  obj <- NormalizeData(obj,scale.factor = scale_factor)
  obj <- FindVariableFeatures(obj,nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj,features = VariableFeatures(obj))
  pca<-obj@reductions$pca@ cell.embeddings
  pca_15=pca[,1:n_pcs]
  # 6) Distances & silhouette
  d   <- dist(pca_15)
  si <- silhouette(labs, d)    # si[,3] = silhouette width per cell
  
  # 7) Mean ASW per cell type (in lvls order)
  means <- tapply(si[, 3], INDEX = labs, FUN = mean, na.rm = TRUE)
  res <- matrix(as.numeric(means), ncol = 1)
  rownames(res) <- lvls         # ✅ this is correct
  colnames(res) <- "ASW"
  
  return(res)
}

