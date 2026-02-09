## Clustering evaluation.   
Here we use PBMC 10K's Seurat CCA (15 dimension), anchor based transfer label as an example, shows the clustering label evaluation steps.  

**Evaluation**

The resulting clustering evaluation, regardlesss of joint clustering and label transfering, are done on rna and atac label.    
For joint clusterng, we directly evaluate the resulting RNA label and ATAC label.   
For transfer label, we read in the been transferred RNA label(all used Seurat RNA label) with the resulting atac label.


**Load in function and libraries**
```r
source("../Functions/Clustering/Func_label_eval.r")
library(Matrix)
library(stringr)

Seurat_RNA_label<-readRDS("seurat_RNA_label.rds")
atac_label<-readRDS("pbmc_sig_15_atac_label.rds")

cluster_eval<-LB_Eval(Seurat_RNA_label,atac_label)
```


