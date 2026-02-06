#!/usr/bin/env Rscript

suppressMessages({
  library(Seurat)
  library(rliger)
  library(stringr)
  library(Matrix)
})

### For example pbmc10k data
ds       <- "pbmc"
gas_name <- "signac"   # signac / cicero / maestro
k        <- 15

basepath <- "PATH/TO/PROCESSED_DATA"   # user should edit
gaspath  <- "PATH/TO/GAS_RESULTS"      # user should edit

# inputs: sparse matrices (genes x cells)
rna <- readRDS(file.path(basepath, ds, "dgCMatrix_count_rna.rds"))
gas <- readRDS(file.path(gaspath, sprintf("%s_%s.rds", ds, gas_name)))

# cicero GAS needs rescaling
if (gas_name == "cicero") {
  gas@x <- gas@x * 1e4
}

# prefix RNA cell IDs to keep modalities distinguishable
colnames(rna) <- paste0("rna_", colnames(rna))

# build LIGER object
rna <- as.matrix(rna)
gas <- as.matrix(gas)

data_list <- list(
  atac = gas,
  rna  = rna
)

int.data <- createLiger(data_list)

rm(data_list, gas, rna)
gc()

# LIGER pipeline
set.seed(1)
int.data <- normalize(int.data)
int.data <- selectGenes(int.data, datasets.use = c("atac", "rna"),
                        nGenes = c(3000, 1000))
int.data <- scaleNotCenter(int.data)
int.data <- runIntegration(int.data, k = k)
int.data <- quantileNorm(int.data)
int.data <- runCluster(int.data, resolution = 0.2)

# extract embeddings (cells x k)
H <- int.data@H.norm

rna_idx  <- grepl("^rna_", rownames(H))
atac_idx <- !rna_idx

rna_em  <- H[rna_idx,  , drop = FALSE]
atac_em <- H[atac_idx, , drop = FALSE]

rownames(rna_em) <- str_replace(rownames(rna_em), "^rna_", "")

# extract labels
labs <- int.data@cellMeta$leiden_cluster

rna_label  <- labs[rna_idx]
names(rna_label) <- str_replace(names(rna_label), "^rna_", "")

atac_label <- labs[atac_idx]

# outputs
utils::write.csv(
  rna_em,
  sprintf("%s_rna_%s_dim%d.csv", ds, gas_name, k),
  quote = FALSE,
  row.names = TRUE
)

utils::write.csv(
  atac_em,
  sprintf("%s_atac_%s_dim%d.csv", ds, gas_name, k),
  quote = FALSE,
  row.names = TRUE
)

label <- list(
  rna_label  = rna_label,
  atac_label = atac_label
)

saveRDS(
  label,
  sprintf("%s_labels_%s_dim%d.rds", ds, gas_name, k)
)
