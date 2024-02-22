## LIGER dimension reduction code   

(pbmc data as an example)  
This script id based on LIGER tutorial :[Joint definition of cell types from single-cell gene expression and chromatin accessibility data](https://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html)   
In Stage I: _Preprocessing and Normalization_,   
1-4 part is calculating gene activity score, which we have made already.  
This script will start from Stage I, step 5 _Normalization_ and the downstream dimension reduction process(stage II), we run with different gene activity scores and choose different numbers of k (size of the reduced dimension)  
To simplify the process, we made a function _Run.liger_ to repeat the dimension reduction steps.
```r
setwd("/data2/duren_lab/naqing/pipeline_building/Dimension_reduction/results/LIGER/")
source("/data2/duren_lab/naqing/pipeline_building/functions/Run.liger.r")
library(Seurat)
library(rliger)
h5_file="/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
counts <- Read10X_h5(h5_file)
pbmc.rna<-counts$'Gene Expression'
# filter cells by barcode
barcode_use<-read.table("/data2/duren_lab/naqing/data/pbmc_10k/barcode_use.txt")
bar<-pbmc.rna@Dimnames[[2]]
idx<-match(barcode_use$V1,bar)
rna<-pbmc.rna[,idx]
gas_names=c("sig","lig","ci","mae")
DIMS=c(15,30,45,60)
gene.activities={}
gene.activities[[1]]<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/pbmc_signac.rds")
gene.activities[[2]]<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/pbmc_liger.rds")
gene.activities[[3]]<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/pbmc_cicero.rds")
gene.activities[[4]]<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/pbmc_maestro.rds")

for (i in 1:4){
    for (j in 1:4){
        embeddings <- Run.liger(rna,gene.activities[[i]],DIMS[j])
        saveRDS(embeddings,paste("pbmc_",gas_names[i],"_",DIMS[j],"embeds.rds"))
      }
    }
