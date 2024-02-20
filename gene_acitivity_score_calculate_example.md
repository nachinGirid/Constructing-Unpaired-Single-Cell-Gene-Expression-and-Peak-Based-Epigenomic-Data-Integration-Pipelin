## Gene_acitivity_score_calculating code, use pbmc10k data as example
### Libraries
```r
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)   #require internet connection
library(cicero)
library(MAESTRO)
library(rliger)
```
### signac
```r
# load the data
counts1 <- Read10X_h5("/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
counts=counts1$Peaks
# annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
# annotation<-readRDS("/data2/duren_lab/naqing/pipeline_building/Data_sets/annotation_files/annotation_EnsDb.Hsapiens.v86.rds") ## pre-made annotation
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments="/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
  min.cells = 0,
  min.features = 0,
)
pbmc.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
  annotation=annotation
)
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 'q0')
pbmc.atac <- RunSVD(pbmc.atac)
gene.activities <- GeneActivity(pbmc.atac)
```
### LIGER
bash
```bash
# go to the folder with fragment files, sor fragment file
cd /data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/
sort -k1,1 -k2,2n -k3,3n pbmc_granulocyte_sorted_10k_atac_fragments.tsv > atac_fragments.sort.bed
# get rid of special chromosomes
awk '$1 !~ /JH|GL|KI/' atac_fragments.sort.bed > filtered_atac_fragments.sort.bed
# sort annotation files
sort -k 1,1 -k2,2n -k3,3n /data2/duren_lab/naqing/pipeline_building/Data_sets/annotation_files/Gene_Annotations_Hsapiens_v86.bed > EnsDb_Hsapiens_v86_genes.sort.bed
sort -k 1,1 -k2,2n -k3,3n /data2/duren_lab/naqing/pipeline_building/Data_sets/annotation_files/Promoter_Annotations_Hsapiens_v86.bed > EnsDb_Hsapiens_v86_promoters.sort.bed
# bedmap command to calculate overlapping elements between indexes and fragment output files:
bedmap --ec --delim '\t' --echo --echo-map-id EnsDb_Hsapiens_v86_promoters.sort.bed filtered_atac_fragments.sort.bed > atac_promoters_bc.bed
bedmap --ec --delim '\t' --echo --echo-map-id EnsDb_Hsapiens_v86_genes.sort.bed filtered_atac_fragments.sort.bed > atac_genes_bc.bed
")
```
R
```r
genes.bc <- read.table(file = "atac_genes_bc.bed", sep = "\t", as.is = c(4,7), header = FALSE)
promoters.bc <- read.table(file = "atac_promoters_bc.bed", sep = "\t", as.is = c(4,7), header = FALSE)
bc <- genes.bc[,7]
bc_split <- strsplit(bc,";")
bc_split_vec <- unlist(bc_split)
bc_unique <- unique(bc_split_vec)
bc_counts <- table(bc_split_vec)
bc_filt <- names(bc_counts)[bc_counts > 1500]
barcodes <- bc_filt
gene.counts <- makeFeatureMatrix(genes.bc, barcodes)
promoter.counts <- makeFeatureMatrix(promoters.bc, barcodes)
if(dim(promoter.counts)[1]>dim(gene.counts)[1]){
promoter.counts <- promoter.counts[rownames(gene.counts),]
}
if(dim(promoter.counts)[1]<dim(gene.counts)[1]){
gene.counts <- gene.counts[rownames(promoter.counts),]
}
atac <- gene.counts + promoter.counts
liger<-as(atac, "dgCMatrix")

```
### cicero
```r
library(Seurat)
library(monocle3)
library(cicero)
library(Signac)
library(stringr)

# read in matrix data using the Matrix package
"/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
barcode_use<-read.table("/data2/duren_lab/naqing/data/pbmc_10k/barcode_use.txt")
indata<-inputdata.10x$Peaks
bar<-indata@Dimnames[[2]]
idx<-match(barcode_use$V1,bar)
indata<-indata[,idx]
# binarize the matrix
indata@x[indata@x > 0] <- 1

# peankinfo (rows) and cellinfo (columns) used to make a cds object_crucial to calculate GAS
# 1 cellinfo
cellinfo <- data.frame(colnames(indata))
colnames(cellinfo)='V1'
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

#2 peak info
# format peak info
peakinfo=indata@Dimnames[[1]]
peakinfo1=unlist(str_split(peakinfo,':'))
chr0=peakinfo1[(1:length(peakinfo))*2-1]
peakinfo3=unlist(str_split(peakinfo1[(1:length(peakinfo))*2],'-'))
start0=peakinfo3[(1:length(peakinfo))*2-1]
end0=peakinfo3[(1:length(peakinfo))*2]
peakinfo=cbind(chr0,start0,end0)
peakinfo=data.frame(peakinfo)
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name
row.names(indata) <- row.names(peakinfo)

# make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,cell_metadata = cellinfo, gene_metadata = peakinfo))

#saved 
#input_cds <- readRDS("/data/duren_lab/naqing/Methods_Benchmark/Cicero/input_cds.rds")
input_cds <- monocle3::detect_genes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] #Ensure there are no peaks included with zero reads

#gene_anno
gene_anno <- rtracklayer::readGFF("/data2/duren_lab/naqing/pipeline_building/Data_sets/annotation_files/Homo_sapiens.GRCh38.105.gtf.gz")
gene_anno$chromosome <- gene_anno$seqid
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

### Add a column for the pData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of 
# the genes. We use this as a proxy for promoters.
# To do this we need the first exon of each transcript
pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
# remove all but the first exons per transcript
pos <- pos[!duplicated(pos$transcript),] 
# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1 

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$transcript),] 
neg$start <- neg$end - 1

gene_annotation<- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation[,c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"
gene_annotation_sub <- subset(gene_annotation_sub, !grepl("GL|KI|MT",gene_annotation_sub[,1]))
gene_annotation_sub [,1]<-paste0("chr",gene_annotation_sub[,1] )

input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)
# make conn
hg38<-readRDS("/data2/duren_lab/naqing/pipeline_building/Data_sets/annotation_files/hg38.rds")
conns <- run_cicero(input_cds,hg38) 

#### Generate gene activity scores ####
# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
# remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]
# make a list of num_genes_expressed
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))
# normalize
cicero_gene_activities <- normalize_gene_activities(unnorm_ga,num_genes)
```
### Maestro
```r
library(Seurat)
library(MAESTRO)
library(reticulate)
inputdata.10x <- Read10X_h5("/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
atac<- inputdata.10x$Peaks
barcode_use<-read.table("/data2/duren_lab/naqing/data/pbmc_10k/barcode_use.txt")
bar<-atac@Dimnames[[2]]
idx<-match(barcode_use$V1,bar)
atac<-atac[,idx]
use_python("/data2/duren_lab/naqing/conda_envs/R43/bin/python", required = TRUE)
gas <- ATACCalculateGenescore(atac, organism = "GRCh38")
```
