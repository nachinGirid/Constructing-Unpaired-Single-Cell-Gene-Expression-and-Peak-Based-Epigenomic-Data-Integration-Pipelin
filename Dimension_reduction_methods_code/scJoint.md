## scJoint dimension reduction code  
**1. Prepare input files**  
To Run scJoint main step is straightforward, after the following files are prepared.  
- RNA & ATAC data in .h5 file
- RNA label in a .csv file

**For .h5 files:** Since scJoint provides function _data_to_h5.R_ to prepare a .h5 file from the SingleCellExperiment(sce) object, we can first make sce from our data  

Note, it only use gene activity score as atac data
```r
library(SingleCellExperiment)
library(Seurat)
h5_file="/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
counts <- Read10X_h5(h5_file)
#RNA
pbmc.rna<-counts$'Gene Expression'
barcode_use<-read.table("/data2/duren_lab/naqing/data/pbmc_10k/barcode_use.txt")
bar<-pbmc.rna@Dimnames[[2]]
idx<-match(barcode_use$V1,bar)
pbmc.rna<-pbmc.rna[,idx]

#ATAC
pbmc.atac<- readRDS("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/pbmc_signac.rds")

#create sce
rna<- CreateSeuratObject(counts = pbmc.rna)
atac<- CreateSeuratObject(counts = pbmc.atac)
RNA<-as.SingleCellExperiment(rna)
ATAC<-as.SingleCellExperiment(atac)
```
Now save it as .h5 files.  
```
source("data_to_h5.R")
common_genes <- intersect(rownames(ATAC),
                          rownames(RNA))
length(common_genes)
RNA <- logcounts(RNA[common_genes, ])
ATAC <- logcounts(ATAC[common_genes, ])
write_h5_scJoint(exprs_list = list(rna = RNA,
                                   atac = ATAC), 
                 h5file_list = c("pbmc_rna.h5", 
                                 "pbmc_atac.h5")
                )
```
**For the RNA label in csv file** just read in and check the order, and save as csv.  
```
celltype<-"/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/seurat_object_label/seurat_RNA_label.rds"
barcode_use<-read.table("/data/duren_lab/naqing/data/pbmc_10k/barcode_use.txt")
sum (rownames(celltype)==barcode_use$V1)
sum (rownames(celltype)==colnames(ATAC))
sum (rownames(celltype)==colnames(RNA))
write_csv_scJoint(cellType_list =  list(rna = celltype$x),
                  csv_list = c("rna_label.csv")
                  )
```
