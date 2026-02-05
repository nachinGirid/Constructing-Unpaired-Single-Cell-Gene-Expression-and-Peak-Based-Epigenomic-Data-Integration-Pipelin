# Three-step benchmarking scRNA peak-based Epigenomic Data integration methods
___________________________________________________________________________________________________________________________________________________________________________
## Overview
### Benchmarked done at 3 steps
- Gene activity score method 
- Dimension reduction method
- Clustering and labeling methods
### 12 sets of paired multi-omic data are used in the benchmarking step.  
| Tissue / Dataset | Species | Genome | Cells | Data Modality |
|---|---|---|---:|---|
| PBMC | Human | hg38 | 9,543 | scRNA-seq + scATAC-seq |
| Human brain | Human | hg38 | 3,332 | scRNA-seq + scATAC-seq |
| Mouse brain E18 | Mouse | mm10 | 4,881 | scRNA-seq + scATAC-seq |
| Bone marrow (NeurIPS 2021) | Human | hg38 | 69,249 | scRNA-seq + scATAC-seq |
| Human kidney cancer | Human | hg38 | 22,772 | scRNA-seq + scATAC-seq |
| Small intestine | Human | hg38 | 10,640 | scRNA-seq + scATAC-seq |
| Mouse kidney cancer | Mouse | mm10 | 14,652 | scRNA-seq + scATAC-seq |
| Mouse skin | Mouse | mm10 | 34,774 | scRNA-seq + scATAC-seq |
| Mouse brain | Mouse | mm10 | 3,293 | scRNA-seq + scATAC-seq |
| Mouse frontal cortex & hippocampus (H3K4me1) | Mouse | mm10 | 12,962 | scRNA-seq + histone ChIP-seq |
| Mouse frontal cortex & hippocampus (H3K4me3) | Mouse | mm10 | 7,465 | scRNA-seq + histone ChIP-seq |
| Mouse frontal cortex & hippocampus (H3K27ac) | Mouse | mm10 | 11,749 | scRNA-seq + histone ChIP-seq |


### At each stage, conduct all possible approaches, and evaluate the results using evaluation matrices(Methods).
___________________________________________________________________________________________________________________________________________________________________________
## 1. Gene activity score
#### Methods for Gene Activity Score (GAS)
We have benchmarked 4 gene activity calculating methods each representing different approaches to predicting gene activity from scATAC-seq data. We have followed the online tutorial for each of the methods :
1.  Signac(1.10.0), [Tutorial](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)
2.  LIGER (1.0.0) , [Tutorial](https://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html)
3.  Cicero(1.3.9)   [Tutorial](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#cicero-gene-activity-scores)
4.  MAESTRO(1.5.1) [directly followed the description of the function](https://github.com/liulab-dfci/MAESTRO/blob/master/R/ATACCalculateGenescore.R)
#### Evaluation metrics for gene activity score
1. **Pearson correlation coefficient (PCC)** of each gene's scRNA-seq data and GAS data.
2. **Local neighborhood consistency (LNC)** of each cell's scRNA-seq profile and GAS profile.
3. **Cell-type Average Silhouette width (ASW)** of each cell type on PCA reduced embedding of GAS profile.
___________________________________________________________________________________________________________________________________________________________________________
## 2.Dimension reduction
#### Methods for Dimension Reduction
We benchmarked multiple dimension-reduction methods that embed scRNA-seq and scATAC-seq cells into a shared latent space. For each method, we followed the official tutorial:

1. **Seurat (v4.3.0)** — [Tutorial](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)  
2. **bindSC (v1.0.0)** — [Tutorial](https://htmlpreview.github.io/?https://github.com/KChenlab/bindSC/blob/master/vignettes/mouse_retina/retina.html)  
3. **CoupledNMF (MATLAB implementation)** — [Code & docs](https://web.stanford.edu/group/wonglab/zduren/CoupledNMF/index.html)  
4. **LIGER (v2.1.0)** — [Tutorial](https://welch-lab.github.io/liger/articles/Integrating_scRNA_and_scATAC_data.html)  
5. **scJoint** — [Tutorial](https://github.com/sydneybiox/scJoint/blob/main/tutorial/Analysis%20of%2010xGenomics%20data%20using%20scJoint.ipynb)  
6. **scDART** — [Repo](https://github.com/PeterZZQ/scDART)  
7. **uniPort (v1.2.2)** — [Tutorial](https://uniport.readthedocs.io/en/latest/)  
8. **GLUE (v0.3.2)** — [Tutorial](https://scglue.readthedocs.io/en/latest/tutorials.html)  
9. **MultiMAP** — [Tutorial](https://github.com/Teichlab/MultiMAP)  
10. **SIMBA (v1.2)** — [Tutorial](https://simba-bio.readthedocs.io/en/latest/multiome_10xpmbc10k_integration.html)  
  
#### Evaluation metrics for common embeddings  
1. **Percentage of cells find themselves in 100 nearest neighbors (PFS100NN)** across modalities, take an average of measurements of RNA to ATAC and ATAC to RNA.  
2. **Contrastive similarity** of a cell to its paired cell against the mean similarity to all unpaired cells.  
3. **Cell-type Average Silhouette width (ASW)** of each cell type RNA and ATAC embedding (pbmc data), check both for the ASW on methods' default-size embedding and PCA reduced default-size embeddings to uniformed-sized embeddings.  
___________________________________________________________________________________________________________________________________________________________________________
## 3.Clustering and labeling 
#### Methods for Clustering and labeling  
We benchmarked multiple strategies for clustering and transferring RNA-derived labels to ATAC cells:

1. **Seurat anchor-based label transfer** — using `FindTransferAnchors` and `TransferData` with default parameters  
2. **FigR anchor-based label transfer** — [Tutorial](https://buenrostrolab.github.io/FigR/articles/FigR_stim.html)  
3. **Leiden clustering** — performed in Scanpy on the joint embedding (cosine distance; neighbors = 20, except BMMC = 200; resolution = 1.0)  
4. **KNN label transfer** — `KNeighborsClassifier` (Euclidean distance; k = 5)  
5. **Optimal transport (OT) label transfer** — Moscot (v0.4.0) — [Tutorial](https://moscot.readthedocs.io/en/latest/notebooks/tutorials/600_tutorial_translation.html#)


#### Evaluation metrics for clustering  

1. **Normalized mutual information(NMI):** between RNA label and ATAC label  
2. **Average label consistency** percentage of label match, per cluster ( average of two percentages: for each cluster RNA cluster, percentage of matched ATAC cells with the same label. Calculate the same for ATAC. Then, we get the mean of all RNA clusters and all ATAC clusters. )
___________________________________________________________________________________________________________________________________________________________________________

# Script examples (PBMC 10k as an example dataset)

This section provides end-to-end example scripts using **PBMC 10k** to demonstrate:
- Running each method (GAS → dimension reduction → clustering/labeling)
- Evaluating outputs

___________________________________________________________________________________________________________________________________________________________________________

## Running the methods

1. **Gene activity score (GAS) calculation**  
   - [Example script](https://github.com/nachinGirid/Benchmark_scRNA_scATAC_integration_methods/blob/main/gene_acitivity_score_calculate_example.md)

2. **Dimension reduction (joint embedding)**  
   - [Example script (placeholder)](PLACEHOLDER_DIMENSION_REDUCTION_SCRIPT_LINK)  
   - *Planned coverage:* Seurat / LIGER / bindSC / scJoint / scDART / uniPort / GLUE / MultiMAP / SIMBA / CoupledNMF

3. **Clustering & labeling / label transfer**  
   - [Example script (placeholder)](PLACEHOLDER_CLUSTERING_LABELING_SCRIPT_LINK)  
   - *Planned coverage:* Seurat anchors / FigR / Leiden / KNN / OT (Moscot)

___________________________________________________________________________________________________________________________________________________________________________

## Evaluation

1. **GAS evaluation**  
   - [Script for GAS evaluation](https://github.com/nachinGirid/Benchmark_scRNA_scATAC_integration_methods/blob/main/gene_activity_score_evaluation.md)

2. **Integration evaluation (placeholder)**  
   - [Example script (placeholder)](PLACEHOLDER_INTEGRATION_EVALUATION_SCRIPT_LINK)

3. **Clustering / labeling evaluation (placeholder)**  
   - [Example script (placeholder)](PLACEHOLDER_CLUSTERING_EVALUATION_SCRIPT_LINK)

4. **Evaluation functions (metrics & utilities)**  
   - [Evaluation functions (placeholder)](PLACEHOLDER_EVALUATION_FUNCTIONS_LINK)  
   - *Planned content:* reusable metric functions (e.g., paired-cell MNN, label transfer accuracy, mixing metrics), plotting helpers, and config/constants.

___________________________________________________________________________________________________________________________________________________________________________

