# Seurat (RNA labels + preprocessing objects)

This folder contains the **Seurat preprocessing and Seurat CCA scripts** used in our benchmark pipeline.  
They provide two reusable inputs for downstream steps (Seurat object and RNA label that will be used as a reference label to be transferred by the label transfer method in the clustering/labeling step).

We show how to:

1) **Generate Seurat RNA labels** by clustering RNA alone.  

2) **Preprocess RNA and ATAC into Seurat objects** from raw count matrices.  
   These preprocessed Seurat objects are reused as standardized inputs for multiple integration methods implemented in this repository.

---

## Scripts (run in this order)

### 1) Generate RNA labels (RNA-only clustering)

**Script:** `scripts/00_make_seurat_rna_labels.R`

**Input**
- `data/pbmc_rna_raw_counts.rds`  
  Raw RNA count matrix (genes × cells), stored as `dgCMatrix` (recommended).

**Output**
- `processed_data/pbmc/seurat_RNA_label.rds`  
  A named factor of Seurat cluster IDs (names = RNA cell barcodes).

**Run**
```bash
Rscript scripts/00_make_seurat_rna_labels.R \
  data/pbmc_rna_raw_counts.rds \
  processed_data/pbmc/seurat_RNA_label.rds
