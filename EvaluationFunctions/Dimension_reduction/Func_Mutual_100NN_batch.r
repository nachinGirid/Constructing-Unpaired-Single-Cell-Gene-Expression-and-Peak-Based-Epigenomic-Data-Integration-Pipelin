Mutual_100NN_batched <- function(rna_embedding, atac_embedding, method = "cosine", batch_size = 7000) {
  library(proxy)  # For distance calculation
  library(Matrix)
  set.seed(123)  # For reproducibility
  rna_embedding=rna_embedding[rowSums(rna_embedding)!=0,]
  atac_embedding=atac_embedding[rowSums(atac_embedding)!=0,]
  common_cell<-intersect(rownames(rna_embedding),rownames(atac_embedding))
  atac_embedding=atac_embedding[common_cell,]
  rna_embedding=rna_embedding[common_cell,] 
  ncell=dim(rna_embedding)[1]
  k <- round(ncell * 0.01)
  
  stopifnot(sum(rownames(atac_embedding) == rownames(rna_embedding)) == ncell)
  
  # Initialize binary vector for mutual nearest neighbors
  mnn_binary <- rep(0, ncell)

  print("Calculating distances in batches...")
  
  # Process data in batches
  for (start in seq(1, ncell, by = batch_size)) {
    end <- min(start + batch_size - 1, ncell)
    
    # Extract the batch
    rna_batch <- rna_embedding[start:end, , drop = FALSE]
    atac_batch <- atac_embedding[start:end, , drop = FALSE]
    
    # Calculate distances for the current batch
    dis_batch_RA <- proxy::dist(rna_batch, atac_embedding, method = method)
    dis_batch_AR <- proxy::dist(atac_batch, rna_embedding, method = method)
    dis_batch_RA <- as.matrix(dis_batch_RA)
    dis_batch_AR <- as.matrix(dis_batch_AR)
    
    # Calculate nearest neighbors for the batch
    for (i in 1:nrow(dis_batch_RA)) {
      global_idx <- start + i - 1
      
      # Top k neighbors in ATAC for RNA cell
      ORDER_RA <- order(dis_batch_RA[i, ], decreasing = FALSE)[1:k]
      
      # Top k neighbors in RNA for ATAC cell
      ORDER_AR <- order(dis_batch_AR[i, ], decreasing = FALSE)[1:k]
      
      # Record mutual neighbors
      if (global_idx %in% ORDER_RA & global_idx %in% ORDER_AR) {
        mnn_binary[global_idx] <- 1
      }
    }
  }
  
  # Calculate overlap of mutual nearest neighbors
  overlap <- sum(mnn_binary)
  
  print("Batch distance calculation complete.")
  
  # Percentage of mutual nearest neighbors
  per <- overlap / ncell
  return(per)
}
