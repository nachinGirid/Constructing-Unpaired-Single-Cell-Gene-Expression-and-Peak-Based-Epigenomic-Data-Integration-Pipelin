## Seurat dimension reduction code
(pbmc data as example)

This script is based on the Seurat Tutorial [Integrating scRNA-seq and scATAC-seq data](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)  
To simplify the process of using different gene activity scores and different dimensions, we made a function _Run.seurat.DMR_, which takes Seuart preprocessed RNA and ATAC Seurat object and corresponding gene activity score and requested size of embedding in, and output reduced RNA and ATAC embeddings.

-------------------------------------------------------------------------------------------------------------------------------
libraries
```r
library(Seurat)
library(Signac)
packageVersion("Seurat")

h5_file="/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
counts <- Read10X_h5(h5_file)
pbmc.rna<-counts$'Gene Expression'
pbmc.atac<-counts$Peak

# filter cells by barcode
barcode_use<-read.table("/data2/duren_lab/naqing/data/pbmc_10k/barcode_use.txt")

bar<-pbmc.rna@Dimnames[[2]]
idx<-match(barcode_use$V1,bar)
pbmc.rna<-pbmc.rna[,idx]

bar<-pbmc.atac@Dimnames[[2]]
idx<-match(barcode_use$V1,bar)
pbmc.atac<-pbmc.atac[,idx]

# make Seurat obj
rna<-CreateSeuratObject(counts=pbmc.rna)
atac<-CreateSeuratObject(counts=pbmc.atac)

# Perform standard analysis of each modality independently 
# RNA analysis
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- RunUMAP(rna, dims = 1:30)

# We exclude the first dimension as this is typically correlated with sequencing depth
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

saveRDS(rna,"pbmc_Seurat_preprocess_obj.rds")
saveRDS(atac,"pbmc_Seurat_preprocess_obj.rds")
```
Once we have preprocessed rna and atac objects, we can repeate following steps to conduct dimension reduction with different gene activity score and with different numbers of dimensions
```
gas_names=c("sig","lig","ci","mae")
DIMS=c(15,30,45,60)
gene.activities={}
gene.activities[[1]]<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/pbmc_signac.rds")
gene.activities[[2]]<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/pbmc_liger.rds")
gene.activities[[3]]<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/pbmc_cicero.rds")
gene.activities[[4]]<-readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/pbmc_maestro.rds")

for (i in 1:4){
    for (j in 1:4){
        embeddings <- Run.seurat.DMR(rna,atac,gene.activities[[i]],1:DIMS[j])
        saveRDS(embeddings,paste("pbmc_",gas_names[i],"_",DIMs[j],"embeds.rds"))
      }
    }
```
