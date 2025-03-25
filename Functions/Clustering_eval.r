# Function to find indices
count_matched <- function(mat) {
  # Initialize a vector to store the indices
  mathced <- numeric(length = nrow(mat))
  
  # Iterate over each row
  for (i in 1:nrow(mat)) {
    row_name <- rownames(mat)[i]
    # Check if the row name exists in the column names
    if(row_name %in% colnames(mat)) {
      # If yes, find the index of the column
      mathced[i] <- mat[i,which(colnames(mat) == row_name)]
    } else {
      # If no match is found, assign 0
      mathced[i] <- 0
    }
  }
  return(mathced)
}

LB_Eval<- function(rna_lb,atac_lb,true_label=NULL){
# IMPORTANT NOTE: rna_lb and atac_lb MUST be numeric and start from 1
    if (!is.numeric(rna_lb)) {
    stop("error: rna_lb must be numeric.")
    }
    if (!is.numeric(atac_lb)) {
    stop("error: atac_lb must be numeric.")
    }
    
    library(aricode)
    library(cluster)
    library(dynutils)
    library(fossil)

    N=length(rna_lb) #cell number

    # NMI,ARI,consist, avrg-consistancy
    label_evaluation=matrix(0,1,4)
    colnames(label_evaluation)<-c('NMI','ARI','over_all_consistancy','Avrg_consistancy')
    label_evaluation[1,1]<-NMI(rna_lb,atac_lb)
    label_evaluation[1,2]<-rand.index(rna_lb,atac_lb)
    label_evaluation[1,3]<-sum(rna_lb==atac_lb)/N
    # cluster label matcing
    rna_clusters<-sort(unique(rna_lb))
    atac_clusters<-sort(unique(atac_lb))    
    tab_rna_atac=matrix(0,length(rna_clusters),length(atac_clusters))
    rownames(tab_rna_atac)<-paste0("Cluster",rna_clusters)
    colnames(tab_rna_atac)<-paste0("Cluster",atac_clusters)
    for (i in 1:length(rna_clusters)){
        for (j in 1:length(atac_clusters)){
         tab_rna_atac[i,j]<-(sum(atac_lb[rna_lb==rna_clusters[i]]==atac_clusters[j]))/(sum(rna_lb==rna_clusters[i]))
        }
    }
    # when atac label
    tab_atac_rna=matrix(0,length(atac_clusters),length(rna_clusters))
    rownames(tab_atac_rna)<-paste0("Cluster",atac_clusters)
    colnames(tab_atac_rna)<-paste0("Cluster",rna_clusters)
    for (i in 1:length(atac_clusters)){
        for (j in 1:length(rna_clusters)){
         tab_atac_rna[i,j]<-(sum(rna_lb[atac_lb==atac_clusters[i]]==rna_clusters[j]))/(sum(atac_lb==atac_clusters[i]))
        }
    }
    # Avrg
    rna_match=mean(count_matched(tab_rna_atac))
    atac_match=mean(count_matched(tab_atac_rna))
    label_evaluation[1,4]<-mean(rna_match,atac_match)   
    return(label_evaluation)
}
