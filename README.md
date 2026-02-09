# Three-step benchmarking scRNA peak-based Epigenomic Data integration methods
___________________________________________________________________________________________________________________________________________________________________________
## Overview
### Benchmarking is done in 3 steps
- Gene activity score method 
- Dimension reduction method
- Clustering and labeling methods
### At each stage, conduct all possible combinations from different steps and evaluate the results using evaluation matrices(Methods).
___________________________________________________________________________________________________________________________________________________________________________
## 1. Gene activity score
#### Methods for Gene Activity Score (GAS)
We have benchmarked 4 gene activity calculating methods each representing different approaches to predicting gene activity from scATAC-seq data. We have followed the online tutorial for each of the methods :
1.  Signac(1.10.0), [Tutorial](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)
2.  LIGER (v2.1.0), [Tutorial](https://welch-lab.github.io/liger/articles/Integrating_scRNA_and_scATAC_data.html)  
3.  Cicero(1.3.9) ,  [Tutorial](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#cicero-gene-activity-scores)
4.  MAESTRO(1.5.1), [directly followed the description of the function](https://github.com/liulab-dfci/MAESTRO/blob/master/R/ATACCalculateGenescore.R)
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
1. **Percentage of same cell's different profiles are mutually in the top 1% nearest neighbors fro each other(MNN%)**.  
2. **Contrastive similarity** of a cell to its paired cell against the mean similarity to all unpaired cells.  
3. **Cell-type Average Silhouette width (ASW)** of each cell type in RNA and ATAC embedding. 
___________________________________________________________________________________________________________________________________________________________________________
## 3.Clustering and labeling 
#### Methods for Clustering and Labeling  
We benchmarked multiple strategies for clustering and transfer label:

1. **Seurat anchor-based label transfer** — using `FindTransferAnchors` and `TransferData` with default parameters  
2. **FigR anchor-based label transfer** — [Tutorial](https://buenrostrolab.github.io/FigR/articles/FigR_stim.html)  
3. **Leiden clustering** — performed in Scanpy on the joint embedding (cosine distance; neighbors = 20, except BMMC = 200; resolution = 1.0)  
4. **KNN label transfer** — `KNeighborsClassifier` (Euclidean distance; k = 5)  
5. **Optimal transport (OT) label transfer** — Moscot (v0.4.0) — [Tutorial](https://moscot.readthedocs.io/en/latest/notebooks/tutorials/600_tutorial_translation.html#)


#### Evaluation metrics for clustering  

1. **Normalized mutual information(NMI):** between RNA label and ATAC label  
2. **Average label consistency** percentage of label match, per cluster ( average of two percentages: for each cluster RNA cluster, percentage of matched ATAC cells with the same label. Calculate the same for ATAC. Then, we get the mean of all RNA clusters and all ATAC clusters. )
___________________________________________________________________________________________________________________________________________________________________________

# Script examples 

This section provides example scripts to demonstrate:
- Running each method (GAS → dimension reduction → clustering/labeling)
- Evaluating outputs

___________________________________________________________________________________________________________________________________________________________________________

## Running the methods

1. **Gene activity score (GAS) calculation**  
   - [Example script](https://github.com/nachinGirid/Constructing-Unpaired-Single-Cell-Gene-Expression-and-Peak-Based-Epigenomic-Data-Integration-Pipelin/tree/main/Run_method/Run_Gene_activity_score)

2. **Dimension reduction (joint embedding)**  
   - [Example script](https://github.com/nachinGirid/Constructing-Unpaired-Single-Cell-Gene-Expression-and-Peak-Based-Epigenomic-Data-Integration-Pipelin/tree/main/Run_method/Run_Dimension_reduction)  
   - *coverage:* Seurat / LIGER / bindSC / scJoint / scDART / uniPort / GLUE / MultiMAP / SIMBA / CoupledNMF

3. **Clustering & labeling / label transfer**  
   - [Example script](https://github.com/nachinGirid/Constructing-Unpaired-Single-Cell-Gene-Expression-and-Peak-Based-Epigenomic-Data-Integration-Pipelin/tree/main/Run_method/Run_Clustering)
   - *coverage:*   OT (Moscot)/FigR /Leiden / KNN / 

___________________________________________________________________________________________________________________________________________________________________________

## Evaluation

1. **GAS evaluation**  
   - [Script for GAS evaluation](https://github.com/nachinGirid/Constructing-Unpaired-Single-Cell-Gene-Expression-and-Peak-Based-Epigenomic-Data-Integration-Pipelin/tree/main/Evaluation/Example_GAS_evaluation)

2. **Integration evaluation**  
   - [Example script](https://github.com/nachinGirid/Constructing-Unpaired-Single-Cell-Gene-Expression-and-Peak-Based-Epigenomic-Data-Integration-Pipelin/tree/main/Evaluation/Example_dimensioin_reduction_evaluation)

3. **Clustering/labeling evaluation**  
   - [Example script](https://github.com/nachinGirid/Constructing-Unpaired-Single-Cell-Gene-Expression-and-Peak-Based-Epigenomic-Data-Integration-Pipelin/tree/main/Evaluation/Example_clustering_evaluation)

4. **Evaluation functions (metrics & utilities)**  
   - [Evaluation functions](https://github.com/nachinGirid/Constructing-Unpaired-Single-Cell-Gene-Expression-and-Peak-Based-Epigenomic-Data-Integration-Pipelin/tree/main/Evaluation/Functions)
___________________________________________________________________________________________________________________________________________________________________________

