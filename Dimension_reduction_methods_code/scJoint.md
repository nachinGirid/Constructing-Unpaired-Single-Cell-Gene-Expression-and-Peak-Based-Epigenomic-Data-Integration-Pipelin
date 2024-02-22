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
**2. Generate .npz files.** This step is done on python
```python
import process_db
import h5py
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import random
random.seed(1)
rna_h5_files = ["data_10x/pbmc_rna.h5"] 
rna_label_files = ["data_10x/rna_label.csv"] # csv file

atac_h5_files = ["data_10x/pbmc_atac.h5"]
atac_label_files = []

process_db.data_parsing(rna_h5_files, atac_h5_files)
rna_label = pd.read_csv(rna_label_files[0], index_col = 0)
rna_label
print(rna_label.value_counts(sort = False))
process_db.label_parsing(rna_label_files, atac_label_files)
```
**3. Edit config.py** Before main.py     
Give the following information:  
1.  Number of clusters
2.  Number of common genes
3.  The paths to the input file
```py
self.number_of_class = 11
self.input_size = 15463
self.rna_paths = ['data_10x/pbmc_rna.npz']
self.rna_labels = ['data_10x/rna_label.txt']		
self.atac_paths = ['data_10x/pbmc_atac.npz']
```

**4. Run main.py** 
