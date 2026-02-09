transfer_rna_labels_to_atac <- function(pairing, rna_lab, atac_pcs, neighbors_number = 5) {
  stopifnot(all(c("ATAC","RNA") %in% colnames(pairing)))
  stopifnot(!is.null(rownames(atac_pcs)))
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Install 'FNN'.")

  # keep names
  if (is.factor(rna_lab)) {
    rna_lab <- setNames(as.character(rna_lab), names(rna_lab))
  }
    rna_lab <- setNames(as.vector(rna_lab), names(rna_lab))
  print(paste("the transfered RNA_label be like ..",head(rna_lab)))
  # 1) direct transfer: RNA -> ATAC
  library(stringr)
  rownames(atac_pcs)<-str_replace(rownames(atac_pcs),'_atac','')
  pairing$RNA<-str_replace(pairing$RNA,'_rna','')
  pairing$ATAC<-str_replace(pairing$ATAC,'_atac','')
  rna_bc <- sub("(_rna)?$", "", pairing$RNA)
  labs   <- rna_lab[rna_bc]
  keep   <- !is.na(labs)
  atac_lab <- stats::setNames(as.character(labs[keep]), pairing$ATAC[keep])

  # 2) fill unmatched ATAC by kNN vote among already labeled ATAC
  all_atac <- rownames(atac_pcs)
  atac_lab <- atac_lab[!is.na(match(names(atac_lab), all_atac))]
  um <- setdiff(all_atac, names(atac_lab))
  if (length(um) > 0 && length(atac_lab) > 0) {
    ref_idx <- match(names(atac_lab), all_atac)
    qry_idx <- match(um, all_atac)
    X_ref <- atac_pcs[ref_idx, , drop = FALSE]
    X_qry <- atac_pcs[qry_idx, , drop = FALSE]
    k_use <- max(1L, min(neighbors_number, nrow(X_ref)))
    kn <- FNN::get.knnx(X_ref, X_qry, k = k_use)
    nn_index <- if (k_use == 1) matrix(kn$nn.index, ncol = 1) else kn$nn.index
    nn_dist  <- if (k_use == 1) matrix(kn$nn.dist,  ncol = 1) else kn$nn.dist
    ref_labels <- unname(as.character(atac_lab))
    voted <- vapply(seq_along(um), function(i){
      idx <- nn_index[i,]; d <- nn_dist[i,]; labs <- ref_labels[idx]
      tab <- table(labs); top <- names(tab)[tab == max(tab)]
      if (length(top) == 1) top else {
        sums <- sapply(top, function(L) sum(d[labs == L], na.rm = TRUE))
        top[which.min(sums)]
      }
    }, character(1))
    atac_lab <- c(atac_lab, stats::setNames(voted, um))
  }

  return(atac_lab[all_atac])
}
