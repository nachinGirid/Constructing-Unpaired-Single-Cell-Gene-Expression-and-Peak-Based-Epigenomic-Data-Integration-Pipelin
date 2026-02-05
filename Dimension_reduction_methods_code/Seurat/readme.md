# Seurat (Dimension reduction, including making the preprocessing Seurat objects and Seurat RNA label generating)

This folder contains the **Seurat preprocessing and Seurat CCA + clustering scripts** used in our benchmark pipeline.  
This folder generates two reusable inputs for downstream steps 
- RNA and ATAC Seurat object
- RNA label that will be used as a reference label to be transferred by the label transfer method in the clustering/labeling step.
---

## Scripts (run in this order)

### 1) Generate RNA labels (RNA-only clustering)
`00_make_seurat_rna_labels.R`  
`01_preprocess_seurat_objects.R`  
`02_cca_transfer.R`

**Output**
- seurat_RNA_label.rds
- seurat_obj_RNA.rds
- seurat-obj_ATAC.rds
- CCA_RNA.cvs
- CCA_ATAC.cvs
