LB_Eval <- function(rna_lb, atac_lb) {
  # deps (same as your originals)
  if (!requireNamespace("aricode", quietly = TRUE)) {
    stop("Please install aricode: install.packages('aricode')")
  }

  # same cells; drop NA pairs
  stopifnot(length(rna_lb) == length(atac_lb))
  keep <- !(is.na(rna_lb) | is.na(atac_lb))
  if (!all(keep)) {
    warning(sprintf("Dropping %d cells with NA labels.", sum(!keep)))
    rna_lb  <- rna_lb[keep]
    atac_lb <- atac_lb[keep]
  }
  N <- length(rna_lb)

  # --- your strict overall consistency (exact equality of original labels) ---
  over_all_consistancy <- sum(rna_lb == atac_lb) / N

  # --- NMI / Rand (label-invariant encodings) ---
  rseq <- as.integer(factor(rna_lb))
  aseq <- as.integer(factor(atac_lb))
  NMI <- aricode::NMI(rseq, aseq)
  ARI <- aricode::ARI(rseq, aseq)  

  # --- build conditional tables EXACTLY as in your old function ---
  rna_clusters  <- sort(unique(rna_lb))
  atac_clusters <- sort(unique(atac_lb))

  tab_rna_atac <- matrix(0, length(rna_clusters), length(atac_clusters))
  rownames(tab_rna_atac) <- paste0("Cluster", rna_clusters)
  colnames(tab_rna_atac) <- paste0("Cluster", atac_clusters)
  for (i in seq_along(rna_clusters)) {
    r <- rna_clusters[i]
    in_r <- (rna_lb == r)
    denom <- sum(in_r)
    for (j in seq_along(atac_clusters)) {
      a <- atac_clusters[j]
      tab_rna_atac[i, j] <- if (denom == 0) 0 else sum(atac_lb[in_r] == a) / denom
    }
  }

  tab_atac_rna <- matrix(0, length(atac_clusters), length(rna_clusters))
  rownames(tab_atac_rna) <- paste0("Cluster", atac_clusters)
  colnames(tab_atac_rna) <- paste0("Cluster", rna_clusters)
  for (i in seq_along(atac_clusters)) {
    a <- atac_clusters[i]
    in_a <- (atac_lb == a)
    denom <- sum(in_a)
    for (j in seq_along(rna_clusters)) {
      r <- rna_clusters[j]
      tab_atac_rna[i, j] <- if (denom == 0) 0 else sum(rna_lb[in_a] == r) / denom
    }
  }

  # --- your diagonal-only averaging (vectorized version of count_matched) ---
  count_matched <- function(mat) {
    rn <- rownames(mat); cn <- colnames(mat)
    idx <- match(rn, cn)                # position of matching column for each row (NA if none)
    out <- numeric(length(rn))
    hit <- !is.na(idx)
    out[hit] <- mat[cbind(which(hit), idx[hit])]
    out
  }

  rna_match  <- mean(count_matched(tab_rna_atac))
  atac_match <- mean(count_matched(tab_atac_rna))
  Avrg_consistancy <- mean(c(rna_match, atac_match))  

  data.frame(
    NMI = NMI,
    ARI = ARI,
    over_all_consistancy = over_all_consistancy,
    Avrg_consistancy = Avrg_consistancy,
    check.names = FALSE
  )
}