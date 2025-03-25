contrastive_sim<-function(dim_rna,dim_atac,temperature=1){
    library(Matrix)
    ncell=dim(dim_rna)[1]
    dd_sparse <- sparseMatrix(i = 1:ncell, j = 1:ncell, x = rep(1, ncell), dims = c(ncell, ncell))
    zero_sparse_matrix <- sparseMatrix(i = integer(0), j = integer(0), dims = c(ncell, ncell))
    DD_sparse1=rbind(zero_sparse_matrix,dd_sparse)
    DD_sparse2=rbind(dd_sparse,zero_sparse_matrix)
    labels=cbind(DD_sparse1,DD_sparse2)
    
    #similarity matrix
    combined_matrix<-as.matrix(rbind(dim_rna,dim_atac))
    sim <- combined_matrix / sqrt(rowSums(combined_matrix * combined_matrix))
    sim <- sim %*% t(sim)
    sim=sim/temperature
    rm(combined_matrix,dim_rna,dim_atac)############## remove #######
    gc()
    contrastive_sim_log=c()
    for (i in 1:nrow(sim)){
        SIM=sim[i,-i]
        LAB=labels[i,-i]
        p_sim=exp(SIM[LAB==1])
        n_sim=exp(SIM[LAB==0])
        n_sim_mean=mean(n_sim)
        contrastive_sim_log[i]=log(p_sim/n_sim_mean)
      }
    MCS<-mean(contrastive_sim_log) 
   return(MCS)
    }
