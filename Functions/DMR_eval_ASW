DMR_ASW<- function(rna_embedding,atac_embedding,true_label_num,true_label_str){

  stopifnot(sum(rownames(atac_embedding)==rownames(rna_embedding))==dim(atac_embedding)[1])
  library(stringr)
  library(aricode)
  library(cluster)
  library(dynutils)
  ###cell_type
  true_label <- as.numeric(true_label_num)
  cell_type=unique(true_label_str)
  L=length(cell_type)
  
  celltype_SI=matrix(0,L,2)
  rownames(celltype_SI)<-cell_type
  colnames(celltype_SI)<-c("RNA_embedding_ASW","ATAC_embedding_ASW")
  
  ### SI 
  si_rna=silhouette(true_label,dist(rna_embedding))
  si_atac=silhouette(true_label,dist(atac_embedding))
  
  
  ### cell type average
  for (i in 1:L){
  celltype_SI[i,1]<- mean(si_rna[true_label_str==cell_type[i],3])
  celltype_SI[i,2]<- mean(si_atac[true_label_str==cell_type[i],3])
  }    
return(celltype_SI)
}
