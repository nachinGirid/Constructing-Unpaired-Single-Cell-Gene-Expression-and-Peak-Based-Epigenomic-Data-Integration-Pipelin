## Gene expression presdict from ATAC-seq data
### Libraries
```r
library(Seurat)
library(Signac)
#library(EnsDb.Hsapiens.v86)   #require internet connection
library(cicero)
library(MAESTRO)
library(rliger)
```
### signac
```r
#load the data
counts1 <- Read10X_h5("/data/duren_lab/naqing/data/pbmc_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
# create a Seurat object containing the RNA adata
pbmc.rna <- CreateSeuratObject(
  counts = counts1$`Gene Expression`,
  assay = "RNA"
)
#standard RNA analysis
pbmc.rna<- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)
#load peak information
counts <- Matrix::readMM('matrix_filted_1.mtx')
PeakName=read.table('PeakName.txt',header=FALSE)
barcode=read.table("/data/duren_lab/naqing/data/pbmc_10k/filtered_feature_bc_matrix/barcodes.tsv",header=FALSE)
rownames(counts)=PeakName$V1
colnames(counts)=barcode$V1
### read in premade annotation(code to generate is below)
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotation) <- "UCSC"
# genome(annotation) <- "hg38"
annotation<-readRDS("/data/duren_lab/naqing/Methods_Benchmark/annotations_EnsDb.Hsapiens.v86.RDS")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments="/data/duren_lab/naqing/data/pbmc_10k/filtered_feature_bc_matrix/fragment.tsv.gz",
  min.cells = 0,
  min.features = 0,
)
#ATAC
pbmc.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)
Annotation(pbmc.atac)=annotation
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 'q0')
pbmc.atac <- RunSVD(pbmc.atac)
gene.activities <- GeneActivity(pbmc.atac)
```
*only on pbmc data*
```r
#filter cells
barcode_use<-read.table("/data/duren_lab/naqing/data/pbmc_10k/barcode_use.txt")
barcode<-pbmc.rna@assays$RNA@counts@Dimnames[[2]]
idx<-match(barcode_use$V1,barcode)
pbmc.rna<-pbmc.rna[,idx]
barcode1<-pbmc.atac@assays$peaks@counts@Dimnames[[2]]
idx1<-match(barcode_use$V1,barcode1)
pbmc.atac<-pbmc.atac[,idx1]
# pre_made annotations for hg38
#annotations<-readRDS("/data/duren_lab/naqing/Methods_Benchmark/annotations_EnsDb.Hsapiens.v86.RDS")
load("/data/duren_lab/Kaya/Unpair/application/pbmc/anno_hg38.RData")

#gene.activities <- GeneActivity(pbmc.atac)
```
### LIGER
bash
```bash
sort -k1,1 -k2,2n -k3,3n GSM4138888_scATAC_BMMC_D5T1.fragments.tsv > atac_fragments.sort.bed
sort -k 1,1 -k2,2n -k3,3n hg19_genes.bed > hg19_genes.sort.bed
sort -k 1,1 -k2,2n -k3,3n hg19_promoters.bed > hg19_promoters.sort.bed
```
bedmap
```bedmap
bedmap --ec --delim "\t" --echo --echo-map-id hg19_promoters.sort.bed atac_fragments.sort.bed > atac_promoters_bc.bed
bedmap --ec --delim "\t" --echo --echo-map-id hg19_genes.sort.bed atac_fragments.sort.bed > atac_genes_bc.bed
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
gene.counts <- gene.counts[order(rownames(gene.counts)),]
promoter.counts <- promoter.counts[order(rownames(promoter.counts)),]
D5T1 <- gene.counts + promoter.counts
colnames(D5T1)=paste0("D5T1_",colnames(D5T1))
bmmc.rna <- read10X(sample.dirs = list("/path_to_sample"), sample.names = list("rna"))
```
### cicero
```r
library(Seurat)
library(monocle3)
library(cicero)
library(Signac)
library(stringr)

# read in matrix data using the Matrix package
inputdata.10x <- Read10X_h5("/data2/duren_lab/naqing/data/MouseBrain/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5")
indata<-inputdata.10x$Peaks

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
gene_anno<-readRDS("/data2/duren_lab/naqing/Methods_Benchmark/Cicero/gene_anno_Mus_musculus.GRCm39.110.rds")
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

#### Add a column for the pData table indicating the gene if a peak is a promoter ####
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

gene_annotation_sub <- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"
#gene_annotation_sub<- readRDS("/data/duren_lab/naqing/Methods_Benchmark/Cicero/gene_annotation_sub.rds")
input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)


# make conn
mm10<-readRDS("/data2/duren_lab/naqing/Methods_Benchmark/Cicero/mm10.rds")
conns <- run_cicero(input_cds,mm10) 

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
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
saveRDS(cicero_gene_activities,"/data2/duren_lab/naqing/Benchmark_mouse_brain/GAS/cicero_GAS.rds")
```
### Maestro
```r
inputdata.10x <- Read10X_h5("/data/duren_lab/naqing/data/pbmc_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
atac_counts <- inputdata.10x$Peaks
pbmc.ATAC.RP.res <- ATACCalculateGenescore(inputMat = atac_counts)
```
