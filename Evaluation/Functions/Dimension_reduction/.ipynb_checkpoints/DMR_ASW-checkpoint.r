DMR_ASW <- function(rna_embedding, atac_embedding, true_label) {
  stopifnot(nrow(rna_embedding) == nrow(atac_embedding))
  stopifnot(identical(rownames(rna_embedding), rownames(atac_embedding)))

  suppressPackageStartupMessages({ library(cluster) })

  # labels
  labs  <- as.integer(factor(true_label))   # 1..K
  lvls  <- levels(factor(true_label))       # cell-type names
  if (length(lvls) < 2) stop("Need at least 2 cell types for silhouette.")

  # distances on embeddings (cells in rows)
  d_rna  <- dist(as.matrix(rna_embedding))
  d_atac <- dist(as.matrix(atac_embedding))

  # silhouettes
  si_rna  <- silhouette(labs, d_rna)
  si_atac <- silhouette(labs, d_atac)

  # mean ASW per cell type (aligned to 1..K order)
  m_rna  <- tapply(si_rna[, 3],  INDEX = labs, FUN = mean, na.rm = TRUE)
  m_atac <- tapply(si_atac[, 3], INDEX = labs, FUN = mean, na.rm = TRUE)

  res <- cbind(RNA = as.numeric(m_rna), ATAC = as.numeric(m_atac))
  rownames(res) <- lvls
  res
}
