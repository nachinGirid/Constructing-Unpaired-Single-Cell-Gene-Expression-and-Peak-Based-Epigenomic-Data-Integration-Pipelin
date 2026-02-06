# data-preprocess

This folder contains the **one-time preprocessing** used in our benchmarking project.

In our benchmark, different methods expect different input formats (e.g., raw feature files, sparse count matrices, Seurat objects, or `.h5ad`). For **fairness** and **convenience**, we generate a common set of intermediate inputs **once per dataset** and reuse them across methods. This ensures that downstream comparisons are not affected by method-specific preprocessing differences and avoids repeating the same conversion steps in every method script.

This folder includes:
- data preprocessing code (raw → `dgCMatrix` → Seurat / SCE → `.rds` / `.h5ad`),
- `seurat_rna_label` code (RNA-only clustering labels),
- this README.

---

## Why we do this once

- **Fairness:** all downstream methods start from the same standardized inputs.
- **Convenience:** avoids repeating file conversion and object construction in each method script.
- **Reproducibility:** the saved `.rds` / `.h5ad` intermediates define exactly what each method received.

---

## Outputs produced by preprocessing

Depending on the dataset and raw inputs, preprocessing generates reusable intermediate files such as:

- Raw sparse matrices:
  - `*_rna_raw_dgCMatrix.rds`
  - `*_atac_raw_dgCMatrix.rds`

- Seurat objects:
  - `*_seurat_rna_obj.rds`
  - `*_seurat_atac_obj.rds`

- Filtered Seurat objects (remove all-zero features):
  - `*_Seurat_obj_rna_fil_rowsum0.rds`
  - `*_Seurat_obj_atac_fil_rowsum0.rds`

- SingleCellExperiment objects:
  - `*_SCE_obj_rna_fil_rowsum0.rds`
  - `*_SCE_obj_atac_fil_rowsum0.rds`

- `.h5ad` exports (used by some methods):
  - `rna_fil_rowsum0.h5ad`
  - `atac_fil_rowsum0.h5ad`
  - `GAS.h5ad`

- Preprocessed Seurat objects:
  - `*_Seurat_preprocess_rna_obj.rds`
  - `*_Seurat_preprocess_atac_obj.rds`

- RNA labels (Seurat clustering):
  - `*_seurat_RNA_label.rds`

---

## Preprocessing code used in this project

Below are the scripts used for preprocessing.  
We keep the code in a simple, direct style (no extra folder/path handling beyond what we used in the project).

### A) Create Seurat objects, convert to SCE / H5AD, and preprocess RNA + ATAC

```r
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(zellkonverter)

# Load raw matrices (dgCMatrix)
rna  <- readRDS(rna_path)
atac <- readRDS(atac_path)

# Create Seurat objects
rna  <- CreateSeuratObject(counts = rna)
atac <- CreateSeuratObject(counts = atac)

# Save for downstream use
saveRDS(rna,  "seurat_rna_obj.rds")
saveRDS(atac, "seurat_atac_obj.rds")

# Remove all-zero features
atac_fil <- atac[rowSums(atac) != 0, ]
rna_fil  <- rna[rowSums(rna) != 0, ]

saveRDS(rna_fil,  "Seurat_obj_rna_fil_rowsum0.rds")
saveRDS(atac_fil, "Seurat_obj_atac_fil_rowsum0.rds")

# Convert to SingleCellExperiment (used by some methods)
RNA_sce  <- as.SingleCellExperiment(rna_fil)
ATAC_sce <- as.SingleCellExperiment(atac_fil)

saveRDS(RNA_sce,  "SCE_obj_rna_fil_rowsum0.rds")
saveRDS(ATAC_sce, "SCE_obj_atac_fil_rowsum0.rds")

# Export .h5ad (used by some methods)
zellkonverter::writeH5AD(RNA_sce,  file.path(outdir, "rna_fil_rowsum0.h5ad"))
zellkonverter::writeH5AD(ATAC_sce, file.path(outdir, "atac_fil_rowsum0.h5ad"))

# Standard analysis of each modality independently

## RNA preprocess
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = 5000)
rna <- ScaleData(rna)
rna <- RunPCA(rna)

saveRDS(rna, "Seurat_preprocess_rna_obj.rds")

## ATAC preprocess
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)

saveRDS(atac, "Seurat_preprocess_atac_obj.rds")
```
### B) Generate seurat_rna_label (RNA-only clustering labels)

This script generates RNA labels once per dataset. These labels are reused in all methods that require RNA reference labels (e.g., label transfer).
```r
library(dplyr)
library(Seurat)
library(patchwork)

# Load the raw RNA count as dgCMatrix
rna.data <- readRDS(inpath)

rna <- CreateSeuratObject(counts = rna.data, min.cells = 3)

# Standard preprocess
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)

rna <- RunPCA(rna, features = VariableFeatures(object = rna))

rna <- FindNeighbors(rna, dims = 1:10)
rna <- FindClusters(rna, resolution = 0.5)

saveRDS(rna@active.ident, outpath)

```
