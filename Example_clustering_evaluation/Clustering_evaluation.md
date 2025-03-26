## Clustering evaluation.   
Here we use PBMC 10K's Seurat CCA (15 dimension), anchor based transfer label as an example, shows the clustering label evaluation steps.  

**Evaluation**

The resulting clustering evaluation, regardlesss of joint clustering and label transfering, are done on rna and atac label.    
For joint clusterng, we directly evaluate the resulting RNA label and ATAC label.   
For transfer label, we read in the been transferred RNA label(all used Seurat RNA label) with the resulting atac label.


**Load in function and libraries**
```r
source("Clustering_eval.r")
library(Matrix)
library(stringr)
```
**Read in label**  

*For resluting label from joint clustering method (Leiden)*  

As the Leiden clustering returns a single label for all rna and atac cells, we separate the rna_label and the atac label from the resulting label so taht we can use our evaluation function to evaluate the consistency.
```r
L<-read.table(signac_15_CCA_leiden.txt",sep="\t",header=TRUE)
ncell=length(L$leiden)/2
rna_lb=L$leiden[1:ncell]
atac_lb=L$leiden[(ncell+1 ):length(L$leiden)]
cluster_eval<-LB_Eval(rna_lb,atac_lb)
```
- The *cluster_eval* function returns the NMI and Average consistency together.

 
*For resulting label from label transfering methods*  

We read in the Seurat RNA label, and the atac label that trasferred from it by transfer label aapproaches (Anchor based / KNN / Optimal transport). 
```r
## For example, Seurat anchor based transferred atac label
Seurat_RNA_label<-readRDS("seurat_RNA_label.rds")
RNA_label=as.numeric(Seurat_RNA_label)
atac_labe<-readRDS("pbmc_sig_15_atac_label.rds")
cluster_eval<-LB_Eval(RNA_label,atac_label)
```


