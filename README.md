# Benchmark_scRNA_scATAC_integration_methods
___________________________________________________________________________________________________________________________________________________________________________
## Overview
### Methods are benchmarked in 3 major steps
- Gene activity estimation 
- Dimension reduction
- Clustering and labeling
### Three paired multi-omic data are used as if unpaired to evaluate the performance of the Methods at different stages.  
- pmbc10k
- human brain
- mouse brain
### At each stage, conduct the available approaches, and evaluate the results using certain evaluation matrices.
### Between adjacent stages, test the performance of available combinations of the approaches.
___________________________________________________________________________________________________________________________________________________________________________
## 1. Gene activity score
#### Methods for Gene Activity Score (GAS)
We have benchmarked 4 gene activity calculating methods each representing different approaches to predicting gene activity from scATAC-seq data. We have followed the online tutorial for each of the methods :
1.  _GeneActivity_ Function implemented in Signac(1.10.0), [Tutorial](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)
2.  _makeFeatureMatrix_ Function implemented in LIGER (1.0.0) , [Tutorial](https://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html)
3.  _build_gene_activity_matrix_ Function implemented in Cicero for monocle3 (1.3.9)   [Tutorial](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#cicero-gene-activity-scores)
4.  _ATACCalculateGenescore_ Function implemented in MAESTRO R Package(1.5.1) [directly followed the description of the function](https://github.com/liulab-dfci/MAESTRO/blob/master/R/ATACCalculateGenescore.R)

Here is an [example script](https://github.com/nachinGirid/Benchmark_scRNA_scATAC_integration_methods/blob/main/gene_acitivity_score_calculate_example.md) for making the gene activity scores using these methods. (pbmc data as an example)
#### Evaluation metrics for gene activity score
1. **Pearson correlation coefficient (PCC)** of each gene's scRNA-seq data and GAS data.
2. **Local neighborhood consistency (LNC)** of each cell's scRNA-seq profile and GAS profile.
3. **Cell-type Average Silhouette width (ASW)** of each cell type on PCA reduced embedding of GAS profile.
4. **Normalized correlation** of group-wise gene expression and GAS.

[Script for GAS evaluation](https://github.com/nachinGirid/Benchmark_scRNA_scATAC_integration_methods/blob/main/gene_activity_score_evaluation.md)
___________________________________________________________________________________________________________________________________________________________________________
## 2.Dimension reduction
#### Methods for Dimension Reduction
We have benchmarked. dimension-reduction methods. each method uses different approaches to put scRNA-seq and scATAC-seq cells into a joint latent space.  
We have followed the online tutorial for each of the methods :
1. Seuart() [Tutorial](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)  
2. LIGER()  [Tutorial](https://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html)  
3. scJoint()  [Tutorial](https://github.com/sydneybiox/scJoint/blob/main/tutorial/Analysis%20of%2010xGenomics%20data%20using%20scJoint.ipynb)  
4. Uniport() [Tutorial](https://uniport.readthedocs.io/en/latest/examples/PBMC/pbmc_integration.html)  
5. SIMBA() [Tutorial](https://simba-bio.readthedocs.io/en/latest/multiome_10xpmbc10k_integration.html)  
6. GLUE() [Tutorial](https://scglue.readthedocs.io/en/latest/tutorials.html)  
7. CoupleNMF  
#### Evaluation metrics for common embeddings  
1. **Percentage of cells find themselves in 100 nearest neighbors (PFS100NN)** across modalities, take an average of measurements of RNA to ATAC and ATAC to RNA.  
2. **Contrastive similarity** of a cell to its paired cell against the mean similarity to all unpaired cells.  
3. **Cell-type Average Silhouette width (ASW)** of each cell type RNA and ATAC embedding (pbmc data), check both for the ASW on methods' default-size embedding and PCA reduced default-size embeddings to uniformed-sized embeddings.  
___________________________________________________________________________________________________________________________________________________________________________
## 3.Clustering and labeling 
#### Methods for Clustering and labeling  
1. Seurat  
2. moscot  
3. louvain

#### Evaluation metrics for clustering  
NMI between RNA label and ATAC label
ARI between RNA label and ATAC label
percentage of label match, as whole
percentage of label match, per cluster (calculate for cluster defined by RNA once, and for clusters defined by ATAC again, and get average for each cluster)  

for pbmc data:
NMI between True label with RNA label and ATAC label 

   
 
