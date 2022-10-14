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
# download and unzip
temp <- tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-65/gtf/mus_musculus/Mus_musculus.NCBIM37.65.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match requirements
# input_cds
# read in matrix data using the Matrix package
indata <- Matrix::readMM("filtered_peak_bc_matrix/matrix.mtx") 
# binarize the matrix
indata@x[indata@x > 0] <- 1

# format cell info
cellinfo <- read.table("filtered_peak_bc_matrix/barcodes.tsv")
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

# format peak info
peakinfo <- read.table("filtered_peak_bc_matrix/peaks.bed")
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name
row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
cell_metadata = cellinfo,
gene_metadata = peakinfo))

input_cds <- monocle3::detect_genes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

#gene_anno
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

input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

tail(fData(input_cds))
# DataFrame with 6 rows and 7 columns
#                                 site_name         chr       bp1       bp2
#                                  <factor> <character> <numeric> <numeric>
# chrY_590469_590895     chrY_590469_590895           Y    590469    590895
# chrY_609312_609797     chrY_609312_609797           Y    609312    609797
# chrY_621772_623366     chrY_621772_623366           Y    621772    623366
# chrY_631222_631480     chrY_631222_631480           Y    631222    631480
# chrY_795887_796426     chrY_795887_796426           Y    795887    796426
# chrY_2397419_2397628 chrY_2397419_2397628           Y   2397419   2397628
#                      num_cells_expressed   overlap        gene
#                                <integer> <integer> <character>
# chrY_590469_590895                     5        NA          NA
# chrY_609312_609797                     7        NA          NA
# chrY_621772_623366                   106         2       Ddx3y
# chrY_631222_631480                     2        NA          NA
# chrY_795887_796426                     1         2       Usp9y
# chrY_2397419_2397628                   4        NA          NA

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

# if you had two datasets to normalize, you would pass both:
# num_genes should then include all cells from both sets
unnorm_ga2 <- unnorm_ga
cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), 
                                                    num_genes)
```
### Maestro
```r

inputdata.10x <- Read10X_h5("/data/duren_lab/naqing/data/pbmc_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
atac_counts <- inputdata.10x$Peaks
pbmc.ATAC.RP.res <- ATACCalculateGenescore(inputMat = atac_counts)
```
