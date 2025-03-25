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


PCA_ASW <- function(GAS, true_label) {
  # true_label has to be a vector
  library(Seurat)
  library(qlcMatrix)
  library(aricode)
  library(cluster)

  stopifnot(length(true_label) == dim(GAS)[2])

  # Step 1: Normalize and Scale GAS
  GAS <- log1p(GAS)  # Log-transform (if needed)
  GAS <- scale(GAS)  # Center and scale each feature

  # Step 2: PCA and silhouette calculation
  L = length(unique(true_label))
  pca_SI = matrix(0, L, 1)
  gas_t = t(GAS)
  pca_gas = prcomp(gas_t)
  pc_gas = pca_gas$x
  pc15_gas = pc_gas[, 1:15]
  si = silhouette(true_label, dist(pc15_gas))

  for (i in 1:L) {
    pca_SI[i, 1] <- mean(si[true_label == i, 3])
  }

  return(pca_SI)
}
