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
2.  _makeFeatureMatrix_ Function implemented in rliger (1.0.0) , [Tutorial](https://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html)
3. _build_gene_activity_matrix_ Function implemented in Cicero fro monocle3 (1.3.9)   [Tutorial](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#cicero-gene-activity-scores)
4.  _ATACCalculateGenescore_ Function implemented in MAESTRO R Package(1.5.1) [directly followed the description of the function](https://github.com/liulab-dfci/MAESTRO/blob/master/R/ATACCalculateGenescore.R)

#### Evaluation metrics for gene activity score
1. **Pearson correlation coefficient (PCC)** of each gene's scRNA-seq data and GAS data.
2. **Local neighborhood consistency (LNC)** of each cell's scRNA-seq profile and GAS profile.
3. **Cell-type Average Silhouette width (ASW)** of each cell type on PCA reduced embedding of GAS profile.
___________________________________________________________________________________________________________________________________________________________________________
## 2.Dimension reduction
#### Methods for Dimension Reduction
We have benchmarked. dimension-reduction methods. each method uses different approaches to put scRNA-seq and scATAC-seq cells into a joint latent space.  
We have followed the online tutorial for each of the methods :
Seuart() [Tutorial](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)
LIGER()  [Tutorial](https://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html)
scJoint()  [Tutorial](https://github.com/sydneybiox/scJoint/blob/main/tutorial/Analysis%20of%2010xGenomics%20data%20using%20scJoint.ipynb)
Uniport() [Tutorial]https://uniport.readthedocs.io/en/latest/examples/PBMC/pbmc_integration.html)
SIMBA() [Tutorial](https://simba-bio.readthedocs.io/en/latest/multiome_10xpmbc10k_integration.html)
GLUE() [Tutorial](https://scglue.readthedocs.io/en/latest/tutorials.html)
CoupleNMF

 
