# Dimension-reduction embedding evaluation

**Example:** PBMC 10k — Seurat CCA embedding (15 dimensions) using Signac GAS.  
This README describes steps to prepare and evaluate joint RNA–ATAC embeddings.

---

## Goal
Evaluate the embeddings using three metrics:

1. Mutual nearest-neighbor percentage (MNN)  
2. Average silhouette width (ASW) by cell type  
3. Contrastive similarity

---

## Preprocessing (uniform formatting)

Different tools return embeddings with different shapes and barcode conventions. The preprocessing step normalizes these differences so downstream evaluation is consistent.  



Tasks:
- Remove method-added suffixes from cell barcodes (e.g., `_query`, `_reference`).
- Ensure RNA and ATAC embeddings use identical barcode names and identical row order.
- Save cleaned embeddings to `processed_embeddings/<method>/<dataset>/...` for downstream clustering steps.

**If embeddings are already cleaned and aligned, skip this step.**

Example (Seurat CCA output — run only if needed):

```r
library(Seurat)
em <- readRDS("sample_pbmc_15_cca.rds")

# separate (Seurat stores reference/query rows together)
total_num <- nrow(em)
ncell <- total_num / 2
rna_em  <- em[1:ncell, ]
atac_em <- em[(ncell + 1):total_num, ]

# remove method-added suffixes
library(stringr)
rownames(rna_em)  <- str_replace(rownames(rna_em), "_reference", "")
rownames(atac_em) <- str_replace(rownames(atac_em), "_query", "")

# reorder so RNA and ATAC have identical row order
rna_em <- rna_em[rownames(atac_em), ]
```
## Start evaluation

Source metric functions (located in Functions/Dimension_reduction/):
```r
source("Functions/Dimension_reduction/Func_Mutual_100NN.r")
source("Functions/Dimension_reduction/Func_Contrastive_sim.r")
source("Functions/Dimension_reduction/DMR_ASW.r")
``` 

### Metrics
1) Mutual nearest-neighbor percentage (MNN)

Compute the fraction of cells that are mutual nearest neighbors across modalities:
```R
MNN <- MNN_per(rna_em, atac_em)
```

2) Average silhouette width (ASW) by cell type

Requires ground-truth labels:

```R
sample_label <- readRDS("sample_pbmc_label.rds")
ASW <- DMR_ASW(rna_em, atac_em, sample_label)
```

3) Contrastive similarity

Contrastive measure showing how paired cells are closer than non-paired cells, do not need ground-truth labels:
```r
con_sim <- contrastive_sim(rna_em, atac_em)
```
