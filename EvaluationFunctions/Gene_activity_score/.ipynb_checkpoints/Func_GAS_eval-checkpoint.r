set_common_gene<-function(Gene_expression,GAS_list){
# Gene_expression is gene expression data matrix
# GAS_list is a list contain GAS calculated by different methods
# return a list contain gene expression and different GASs on same set of genes
GAS_common_genes=rownames(GAS_list[[1]])

for(i in 2:length(GAS_list)){
    GAS_common_genes=intersect(GAS_common_genes,rownames(GAS_list[[i]]))
    }
common_genes<-intersect(rownames(Gene_expression),GAS_common_genes)

result=list()
GASs=list()
    
Exp=Gene_expression[common_genes,]
for(i in 1:length(GAS_list)){
    GASs[[i]]<-GAS_list[[i]][common_genes,]
    }
result$Gene_exp<-Exp
result$GASs<-GASs
    return(result)
}

unify_genes_cells <- function(Gene_expression, GAS_list, genes = NULL, cells = NULL) {
  stopifnot(is.matrix(Gene_expression) || inherits(Gene_expression, "Matrix"))
  stopifnot(length(GAS_list) > 0)

  # basic checks
  if (is.null(rownames(Gene_expression)) || is.null(colnames(Gene_expression)))
    stop("Gene_expression must have rownames and colnames.")
  if (any(vapply(GAS_list, function(m) is.null(rownames(m)) || is.null(colnames(m)), TRUE)))
    stop("All GAS matrices must have rownames and colnames.")
  if (any(duplicated(rownames(Gene_expression))) || any(duplicated(colnames(Gene_expression))))
    stop("Duplicate rownames/colnames in Gene_expression are not allowed.")
  if (any(vapply(GAS_list, function(m) any(duplicated(rownames(m))) || any(duplicated(colnames(m))), TRUE)))
    stop("Duplicate rownames/colnames in a GAS matrix are not allowed.")

  # intersect all genes and cells across inputs
  common_genes <- Reduce(intersect, c(list(rownames(Gene_expression)), lapply(GAS_list, rownames)))
  common_cells <- Reduce(intersect, c(list(colnames(Gene_expression)), lapply(GAS_list, colnames)))

  # optional extra filters
  if (!is.null(genes))  common_genes <- intersect(common_genes, genes)
  if (!is.null(cells))  common_cells <- intersect(common_cells, cells)

  # preserve Gene_expression ordering
  common_genes <- rownames(Gene_expression)[rownames(Gene_expression) %in% common_genes]
  common_cells <- colnames(Gene_expression)[colnames(Gene_expression) %in% common_cells]

  # slice
  Gene_expression_aligned <- Gene_expression[common_genes, common_cells, drop = FALSE]
  GAS_aligned <- lapply(GAS_list, function(m) m[common_genes, common_cells, drop = FALSE])
  names(GAS_aligned) <- names(GAS_list)

  list(
    genes = common_genes,
    cells = common_cells,
    Gene_exp = Gene_expression_aligned,
    GASs = GAS_aligned
  )
}

# optional alias
size_unify <- unify_genes_cells



GAS_Exp_Corr<-function(Exp,GAS){

# calculate for one method each time, put gene activity list into one list and run for loop
print("make sure the rows are genes and columns are cells!")
    # libray
    library(Seurat)
    library(qlcMatrix)
    # 
    gene_var<-FindVariableFeatures(Exp)
    var_order<-order(gene_var$vst.variance.standardized,decreasing = TRUE)
    hig_idx<-var_order[1:2000]
    low_idx<-var_order[2001:length(var_order)]
    ### correlation
    a=t(Exp)
    b=t(GAS)
    corr= diag(corSparse(b, Y =a, cov = FALSE))
    stopifnot(length(corr)==dim(Exp)[1])
    corr_high= corr[hig_idx]
    corr_low= corr[low_idx]
    
    corr[is.na(corr)]<-0
    corr_high[is.na(corr_high)]<-0
    corr_low[is.na(corr_low)]<-0

    return(list(corr, corr_high,corr_low))
    }

GAS_Exp_Corr_by_celltype <- function(Exp, GAS, celltype) {
  # calculate for one method each time
  # rows = genes, cols = cells
  message("make sure rows are genes and columns are cells!")

  library(qlcMatrix)

  # basic checks
  stopifnot(ncol(Exp) == ncol(GAS))
  stopifnot(length(celltype) == ncol(Exp))

  # make celltype a factor so we have stable levels
  if (!is.factor(celltype)) {
    celltype <- factor(celltype)
  }
  ct_levels <- levels(celltype)

  n_genes <- nrow(Exp)
  n_ct    <- length(ct_levels)

  gene_names <- rownames(Exp)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(n_genes))
  }

  # output: rows = cell types, cols = genes
  ct_corr <- matrix(NA_real_, nrow = n_ct, ncol = n_genes)
  rownames(ct_corr) <- ct_levels
  colnames(ct_corr) <- gene_names

  for (i in seq_along(ct_levels)) {
    ct <- ct_levels[i]
    idx <- which(celltype == ct)

    # need at least 2 cells to define a correlation
    if (length(idx) < 2) {
      next
    }

    # subset to this cell type
    a <- t(Exp[, idx, drop = FALSE])  # cells × genes (expression)
    b <- t(GAS[, idx, drop = FALSE])  # cells × genes (GAS)

    # pairwise gene–gene correlation; diag is gene-by-gene match
    cmat <- corSparse(b, Y = a, cov = FALSE)
    corr <- diag(cmat)
    stopifnot(length(corr) == n_genes)

    corr[is.na(corr)] <- 0
    ct_corr[i, ] <- corr
  }

  # replace any remaining NAs (e.g. celltypes with <2 cells) with 0
  ct_corr[is.na(ct_corr)] <- 0

  return(ct_corr)
}
                        
                        
                        
                        
                        
                        
LNC<-function(Exp,GAS,block=20){
    #local neighborhood consistancy
    # calculate jaccard similarity to cell itself in two modality
    # could be used to evaluate GAS, how similar cells is when representing by GAS and by gene expression
    # also could be used to evaluate dimension reduction, to see how cross modality similarity is captured for a cell
    # if used for dimension reduction, a transform of matrix might needed
    library(dynutils)
    exp<- as.matrix(Exp)
    gas <- as.matrix(GAS)
    #
    N=dim(exp)[2]
    k=round(N/block)
    #
    RNA_sim<-calculate_similarity(
      x = exp,#RNA
      method = c("cosine"),
      margin = 2,
      diag = FALSE,
      drop0 = FALSE
    )
    ATAC_sim<-calculate_similarity(
      x = gas,# ATAC
      method = c("cosine"),
      margin = 2,
      diag = FALSE,
      drop0 = FALSE
    )
    # minus the expected similarity
    KK=colSums(RNA_sim)
    twom=sum(KK)
    RNA_sim_norm=RNA_sim-(as.matrix(KK) %*% (t(as.matrix(KK/twom))))

    # minus the expected similarity
    KK=colSums(ATAC_sim)
    twom=sum(KK)
    ATAC_sim_norm=ATAC_sim-(as.matrix(KK) %*% (t(as.matrix(KK/twom))))
    

    # KNN for RNA 
    KNN_RNA=matrix(0,N,N)
    for (i in (1:N)) {
        KNN_RNA[,i]=order(RNA_sim_norm[,i],decreasing = TRUE);
      }

    KNN_RNA01=matrix(0,N,N)
    for (i in (1:N)) {
        KNN_RNA01[KNN_RNA[2:k,i],i]=1;
      }

    # KNN for ATAC
    KNN_ATAC=matrix(0,N,N)
    for (i in (1:N)) {
        KNN_ATAC[,i]=order(ATAC_sim_norm[i,],decreasing = TRUE);
      }

    KNN_ATAC01=matrix(0,N,N)
    for (i in (1:N)) {
        KNN_ATAC01[KNN_ATAC[2:k,i],i]=1;
      }
    # shared kNN for both RNA and ATAC
    SNN=t(KNN_RNA01)%*%KNN_ATAC01
    
    self_SNN=c()
    for (i in 1:N){
    self_SNN[i]=SNN[i,i]/k
        }
    
    return(self_SNN)
}



PCA_ASW <- function(GAS, true_label, n_pcs = 15, scale_factor = 1e6) {
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
  
  # 2) Library normalization (CPM-ish, default 1e6)
  csum <- colSums(gas_fil)
  scaling_factor <- csum / scale_factor
  gas_lib_norm <- sweep(gas_fil, 2, scaling_factor, FUN = "/")
  
  # 3) Log-normalization
  gas_norm <- log1p(gas_lib_norm)
  
  # 4) Transpose: rows = cells, cols = genes
  X <- t(gas_norm)
  
  # 5) PCA on cells
  pc <- prcomp(X, center = TRUE, scale. = TRUE)
  k  <- min(ncol(pc$x), n_pcs)
  pcx <- pc$x[, seq_len(k), drop = FALSE]
  
  # 6) Distances & silhouette
  d  <- dist(pcx)
  si <- silhouette(labs, d)    # si[,3] = silhouette width per cell
  
  # 7) Mean ASW per cell type (in lvls order)
  means <- tapply(si[, 3], INDEX = labs, FUN = mean, na.rm = TRUE)
  res <- matrix(as.numeric(means), ncol = 1)
  rownames(res) <- lvls         # ✅ this is correct
  colnames(res) <- "ASW"
  
  return(res)
}


PCA_ASW_cicero <- function(GAS, true_label, n_pcs = 15, log_transform = FALSE) {
  # GAS: genes x cells, *already* normalized by cicero::normalize_gene_activities
  stopifnot(length(true_label) == ncol(GAS))
  suppressPackageStartupMessages(library(cluster))

  # consistent factor mapping
  labs_factor <- factor(true_label)
  labs <- as.integer(labs_factor)   # 1..K
  lvls <- levels(labs_factor)

  # 1) Filter zero-variance genes
  sds <- apply(GAS, 1, stats::sd)
  keep <- !is.na(sds) & sds > 0
  GAS_f <- GAS[keep, , drop = FALSE]

  # 2) Optional log1p (often OK to skip for Cicero)
  if (log_transform) {
    GAS_f <- log1p(GAS_f)
  }

  # 3) Transpose: rows = cells, cols = genes
  X <- t(GAS_f)

  # 4) PCA with centering & scaling
  pc <- prcomp(X, center = TRUE, scale. = TRUE)
  k  <- min(ncol(pc$x), n_pcs)
  pcx <- pc$x[, seq_len(k), drop = FALSE]

  # 5) Distances & silhouette
  d  <- dist(pcx)
  si <- silhouette(labs, d)   # si[,3] = silhouette width per cell

  # 6) Mean ASW per cell type in lvls order
  means <- tapply(si[, 3], INDEX = labs, FUN = mean, na.rm = TRUE)
  res <- matrix(as.numeric(means), ncol = 1)
  rownames(res) <- lvls   # ✅ correct mapping label index -> level name
  colnames(res) <- "ASW"

  return(res)
}

