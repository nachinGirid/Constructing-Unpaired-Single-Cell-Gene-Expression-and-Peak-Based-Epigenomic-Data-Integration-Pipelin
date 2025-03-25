MNN_per<- function(rna_embedding,atac_embedding,method = "cosine"){
  ncell=nrow(atac_embedding)
  k=round(ncell*0.01)
  stopifnot(sum(rownames(atac_embedding)==rownames(rna_embedding))==ncell)
  
  ###  library
  library(proxy)
  library(Matrix)
  ###  Distance
  dis <- as.matrix(proxy::dist(rna_embedding, atac_embedding, method = method))
  N=0
  for(i in 1:ncell){
    ORDER_RA<- order(dis[i,],decreasing = FALSE)[1:k]
    ORDER_AR<- order(dis[,i],decreasing = FALSE)[1:k]
    N=N+((i %in% ORDER_RA)*(i %in% ORDER_AR))
    }
  per=N/ncell
  return(per)
}
