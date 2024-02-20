# Gene activity score evaluation
## Preprocessing
Because different methods generate gene activity scores for different numbers of genes. To compare among four methods, we use the same set of genes. This part of the script reads in gene expression matrix and gene activity scores that are output from different gene activity predicting methods, output a list object contain gene expression and gene activity scores with a common set of genes for downstream evaluation.
**Function and library**
```r
source("/data2/duren_lab/naqing/pipeline_building/functions/Func_GAS_eval.r")
library(Seurat)
library(Matrix)
setwd("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/")
```
**_pbmc 10k_**
```
# Gene expression data
features_matrix<-Read10X_h5("/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Gene_Expression<-features_matrix$'Gene Expression'
barcode_use<-read.table("/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/barcode_use.txt")
barcode<-Gene_Expression@Dimnames[[2]] 
idx<-match(barcode_use$V1,barcode)
Gene_Exp=Gene_Expression[,idx]

# Gene activity scores (GAS)
GASs={}
GASs[[1]]<-readRDS("pbmc_signac.rds")
GASs[[2]]<-readRDS("pbmc_liger.rds")
GASs[[3]]<-readRDS("pbmc_cicero.rds")
GASs[[4]]<-readRDS("pbmc_mastro.rds")

# make a list of gene expression data and gene activity scores with common genes
lst<-set_common_gene(Gene_Exp,GASs)
# save as RDS for later evalution
saveRDS(lst,"/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/common_gene_GASs/pbmc_Exp_5GAS_lst.rds")
```
**_human brain_**
``` R
# Gene expression data
features_matrix<-Read10X_h5("/data2/duren_lab/naqing/pipeline_building/Data_sets/HumanBrain/human_brain_3k_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Gene_Exp<-features_matrix$'Gene Expression'

# Gene activity scores 
GASs={}
GASs[[1]]<-readRDS("humanbrain_signac.rds")
GASs[[2]]<-readRDS("humanbrain_liger.rds")
GASs[[3]]<-readRDS("humanbrain_cicero.rds")
GASs[[4]]<-readRDS("humanbrain_mastro.rds")
# make a list of gene expression data and gene activity scores with common genes
lst<-set_common_gene(Gene_Exp,GASs)
# save as RDS for convinience
saveRDS(lst,"/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/common_gene_GASs/humanbrain_Exp_GAS_lst.rds")
```

**_mouse brain_**
```
# read in gene expression data
features_matrix<-Read10X_h5("/data2/duren_lab/naqing/pipeline_building/Data_sets/MouseBrain/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Gene_Exp<-features_matrix$'Gene Expression'

# read in gene activity scores 
GASs={}
GASs[[1]]<-readRDS("mousebrain_signac.rds")
GASs[[2]]<-readRDS("mousebrain_liger.rds")
GASs[[3]]<-readRDS("mousebrain_cicero.rds")
GASs[[4]]<-readRDS("mousebrain_mastro.rds")
# make a list of gene expression data and gene activity scores with common genes
lst<-set_common_gene(Gene_Exp,GASs)
# save as RDS for convinience
saveRDS(lst,"/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/common_gene_GASs/mousebrain_Exp_GAS_lst.rds")
```
## Evaluation
This part of the script starts to evaluate the gene activity scores with common genes using 3 evaluation matrices.  
**Function and library**
```r
source("/data2/duren_lab/naqing/pipeline_building/functions/Func_GAS_eval.r")
library(Matrix)
setwd("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/evaluation_results/")
```
**Read in results**
```r
#pbmc
lst<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/common_gene_GASs/pbmc_Exp_GAS_lst.rds")
Exp1=lst$Gene_exp
GAS_list1=lst$GASs
#HB
#lst<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/common_gene_GASs/humanbrain_Exp_GAS_lst.rds")
Exp2=lst$Gene_exp
GAS_list2=lst$GASs
#MB
#lst<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/common_gene_GASs/mousebrain_Exp_GAS_lst.rds")
Exp3=lst$Gene_exp
GAS_list3=lst$GASs
```
### 1. Pearson correlation coefficient (PCC) between GAS and gene expression
Empty tables for saving results
```r
# Total 6 tables, 6 box plots. 3 for all genes, 3 for highly variable genes

cor_tab1<-matrix(0,4,dim(Exp1)[1])
cor_tab2<-matrix(0,4,dim(Exp2)[1])
cor_tab3<-matrix(0,4,dim(Exp3)[1])

cor_tab1_hv<-matrix(0,4,2000)
cor_tab2_hv<-matrix(0,4,2000)
cor_tab3_hv<-matrix(0,4,2000)
```
Fill in empty tables with PCC results
```r
for(i in 1:4){
    correlation1<-GAS_Exp_Corr(Exp1,GAS_list1[[i]])
    correlation2<-GAS_Exp_Corr(Exp2,GAS_list2[[i]])
    correlation3<-GAS_Exp_Corr(Exp3,GAS_list3[[i]])

    cor_tab1[i,]<-correlation1[[1]]
    cor_tab2[i,]<-correlation2[[1]]
    cor_tab3[i,]<-correlation3[[1]]

    cor_tab1_hv[i,]<-correlation1[[2]]
    cor_tab2_hv[i,]<-correlation2[[2]]
    cor_tab3_hv[i,]<-correlation3[[2]]
}
# Save evaluation results
corr_tables={}
corr_tables$cor_tab1<-cor_tab1
corr_tables$cor_tab2<-cor_tab2
corr_tables$cor_tab3<-cor_tab3
corr_tables$cor_tab1_hv<-cor_tab1_hv
corr_tables$cor_tab2_hv<-cor_tab2_hv
corr_tables$cor_tab3_hv<-cor_tab3_hv
saveRDS(corr_tables,"corr_tables.rds")
```
### 2. Local neighborhood consistency (LNC) between cells in RNA and GAS profiles
Empty tables for saving results
```r
# Total 3 table, 3 box plots

LNC1<-matrix(0,4,dim(Exp1)[2])
LNC2<-matrix(0,4,dim(Exp2)[2])
LNC3<-matrix(0,4,dim(Exp3)[2])

for(i in 1:4){
    LNC1[i,]<-LNC(Exp1,GAS_list1[[i]],block=20)
    LNC2[i,]<-LNC(Exp2,GAS_list2[[i]],block=20)
    LNC3[i,]<-LNC(Exp3,GAS_list3[[i]],block=20)
}
```
Fill in empty tables with LNC results
```r
# Save evaluation results
LNC_tables={}
LNC_tables$LNC1<-LNC1
LNC_tables$LNC2<-LNC2
LNC_tables$LNC3<-LNC3
saveRDS(LNC_tables,"LNC_tables.rds")
```
### 3. Cell-type Average Silhouette Width (ASW) of PCA reduced GAS profiles (pbmc data)
```r
true_label<-read.table("/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/pbmc_label.txt",header=TRUE,sep="\t")
ASW_GAS=matrix(0,length(unique(true_label$SS)),4)
for (i in 1:4){
ASW_GAS[,i]<-PCA_ASW(Exp,GAS_list[[i]],true_label$SS)
 }
colnames(ASW_GAS)<-c("Signac","liger","cicero","maestro")
saveRDS(ASW_GAS,"ASW_GAS_table.rds")
```
