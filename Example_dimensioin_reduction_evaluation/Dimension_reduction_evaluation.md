## Dimension reduction embedding evaluation.   
Here we use PBMC 10K's Seurat resulting embedding, CCA (15 dimension)using signac GAS, as an example, showing the joint embedding evaluation steps.  

**Preprocessing**  

The resulting dimension reduction  embedding from different methods are in various formats. The purpose of this step is to put embeddings into a uniform, ready-to-evaluate format.    
- Remove the strings added to the cell barcodes by the method during the process of dimension reduction.
- Make the cell embedding barcodes in the same order.
  
As Seurat returns the two modality embeddings in one whole joint embedding matrix, with barcode added "_query " and "_reference" to identify the two modalities, we first separate the two modalities then remove the extra string from barcode.

```r
library(Seurat)
em<-readRDS("sample_pbmc_15_cca.rds")
## separate
total_num=nrow(em)
ncell=total_num/2
rna_em<-em[1:ncell,]
atac_em<-em[(ncell+1):total_num,]
## remove string
library(stringr)
rownames(rna_em)<-str_replace(rownames(rna_em),"_reference","")
rownames(atac_em)<-str_replace(rownames(atac_em),"_query","")
rna_em=rna_em[rownames(atac_em),]
```
**Start evaluation**
```r
source("DMR_eval_MNN_per.r")
source("DMR_eval_ASW.r")
source("DMR_eval_Contrastive_sim.r")
```
### 1. Mutual nearest neighbor percentage  
```
MNN<-MNN_per(rna_em,atac_em)
```
### 2. ASW   
The ground truth cell labels are needed for calculating this metric.  
The numeric label is needed for calculating the silhouette index. 
The string format label is needed when adding the row names for the resulting cell-type average silhouette width table.
```r
sample_label_str<-readRDS("sample_pbmc_str_label.rds")
sample_label_num<-readdRDS("sample_pbmc_num_label.rds")
ASW<-DMR_ASW(rna_em,atac_em,true_label_num,true_label_str)
```
### 3. Contrastive similarity
```r
con_sim<-contrastive_sim(rna_em,atac_em)
```
