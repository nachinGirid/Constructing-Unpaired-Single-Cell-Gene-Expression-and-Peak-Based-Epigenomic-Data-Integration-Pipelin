# GAS (Gene Activity Score) methods

This folder contains code for computing **gene activity scores (GAS)** from scATAC-seq data using multiple commonly used tools:

- **Signac**
- **LIGER**
- **Cicero**
- **MAESTRO**

In our benchmarking project, downstream integration methods expect GAS in a common format (typically a **genes × cells** sparse matrix). We therefore provide method-specific GAS scripts here, while keeping the outputs consistent so they can be used interchangeably in later steps of the pipeline.

**Inputs (typical)**
- `filtered_feature_bc_matrix.h5` (10x multiome feature matrix; optional here)
- `atac_fragments.tsv.gz` (fragment file)
- `atac_raw_dgCMatrix.rds` (ATAC peak-by-cell matrix)
- gene annotation database (`EnsDb.Hsapiens.v86` or `EnsDb.Mmusculus.v79`)

## Signac GAS
```r
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(profmem)

packageVersion("Signac")
packageVersion("Seurat")

h5_file  <- "filtered_feature_bc_matrix.h5"
frag.file <- "atac_fragments.tsv.gz"

# Human example
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- genome

counts <- Read10X_h5(h5_file)

rna  <- readRDS("rna_raw_dgCMatrix.rds")
atac <- readRDS("atac_raw_dgCMatrix.rds")

chrom_assay <- CreateChromatinAssay(
  counts    = atac,
  sep       = c(":", "-"),
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

atac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")

gene.activities <- GeneActivity(atac)
saveRDS(gene.activities, "signac.rds")
```
### LIGER
Required files  
atac_fragments.tsv.gz

hg38_gene_annotation.sort.bed  
hg38_promoter_annotation.sort.bed  
or  
mm10_gene_annotation.sort.bed  
mm10_promoter_annotation.sort.bed    
#### Bash: sort/filter fragments and run bedmap
```bash
atac_fragments.tsv.gz

hg38_gene_annotation.sort.bed
hg38_promoter_annotation.sort.bed
or
mm10_gene_annotation.sort.bed
mm10_promoter_annotation.sort.bed


# --- 2) sort and filter ---
SORTED="tac_fragments.sort.bed"
FILTERED="filtered_atac_fragments.sort.bed"

LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$FRAG_TSV" > "$SORTED"
# drop unplaced/alt contigs like JH*, GL*, KI*
awk '$1 !~ /^(JH|GL|KI)/' "$SORTED" > "$FILTERED"

# --- 3) bedmap overlaps: promoters and genes ---
PROM_OUT="atac_promoters_bc.bed"
GENE_OUT="atac_genes_bc.bed"

# both inputs must be sorted the same way; your annotations already are.
bedmap --ec --delim '\t' --echo --echo-map-id "$PROMS" "$FILTERED" > "$PROM_OUT"
bedmap --ec --delim '\t' --echo --echo-map-id "$GENES" "$FILTERED" > "$GENE_OUT"
```
#### R: build GAS sparse matrix from bedmap outputs  
```r
rna_mat <- reaRDS("rna_raw_dgCMatrix.rds")
genes_bed <- file.path("atac_genes_bc.bed")
proms_bed <- file.path("atac_promoters_bc.bed")
# ---- barcode order from RNA ----
barcode <- colnames(rna_mat)
# ---- read bedmap outputs ----
genes.bc     <- read.table(genes_bed, sep = "\t", header = FALSE, quote = "")
promoters.bc <- read.table(proms_bed, sep = "\t", header = FALSE, quote = "")

# ---- build sparse feature x cell matrices (columns ordered by `barcode`) ----
gene.counts     <- makeFeatureMatrix(genes.bc,     barcode)
promoter.counts <- makeFeatureMatrix(promoters.bc, barcode)

# ---- align features (rows) with union, fill missing with zeros ----
all_feats <- union(rownames(gene.counts), rownames(promoter.counts))

miss_g <- setdiff(all_feats, rownames(gene.counts))
if (length(miss_g)) {
  gene.counts <- rbind(
    gene.counts,
    Matrix(0, nrow = length(miss_g), ncol = ncol(gene.counts), sparse = TRUE,
           dimnames = list(miss_g, colnames(gene.counts)))
  )
}

miss_p <- setdiff(all_feats, rownames(promoter.counts))
if (length(miss_p)) {
  promoter.counts <- rbind(
    promoter.counts,
    Matrix(0, nrow = length(miss_p), ncol = ncol(promoter.counts), sparse = TRUE,
           dimnames = list(miss_p, colnames(promoter.counts)))
  )
}

# enforce identical row order
gene.counts     <- gene.counts[all_feats, , drop = FALSE]
promoter.counts <- promoter.counts[all_feats, , drop = FALSE]

# ---- gene activity ----
atac <- gene.counts + promoter.counts
atac_dgc <- as(atac, "dgCMatrix")

# ---- save ----
saveRDS(atac_dgc,"liger") 
```
### Cicero  
```r
# input
indata <- readRDS("atac_raw_dgCMatrix.rds")

# binarize the matrix
indata@x[indata@x > 0] <- 1

# 1) cellinfo
cellinfo <- data.frame(colnames(indata))
colnames(cellinfo) <- "V1"
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

# 2) peak info
peakinfo <- indata@Dimnames[[1]]
peakinfo1 <- unlist(str_split(peakinfo, ":"))
chr0 <- peakinfo1[(1:length(peakinfo))*2-1]
peakinfo3 <- unlist(str_split(peakinfo1[(1:length(peakinfo))*2], "-"))
start0 <- peakinfo3[(1:length(peakinfo))*2-1]
end0 <- peakinfo3[(1:length(peakinfo))*2]
peakinfo <- cbind(chr0, start0, end0)
peakinfo <- data.frame(peakinfo)
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep = "_")
row.names(peakinfo) <- peakinfo$site_name
row.names(indata) <- row.names(peakinfo)

# make CDS
input_cds <- suppressWarnings(new_cell_data_set(indata, cell_metadata = cellinfo, gene_metadata = peakinfo))
input_cds <- monocle3::detect_genes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0, ]  # keep nonzero peaks

# annotations
if (species == "human") {
  gene_annotation_sub <- readRDS("hg38_gene_annotation.rds")
  genome <- readRDS("human.hg38.genome.rds")
} else {
  gene_annotation_sub <- readRDS("mm10_gene_annotation.rds")
  genome <- readRDS("mouse.mm10.genome.rds")
}

input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

# cicero
conns <- run_cicero(input_cds, genome)
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]

# num_genes_expressed
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

# normalize and save
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
saveRDS(cicero_gene_activities, file.path(("cicero.rds", ds)))
```
### MAESTRO
```r
### maestro##
suppressMessages({
  library(Seurat)
  library(MAESTRO)
  library(reticulate)
})

# python for MAESTRO
use_python("~/.conda/envs/maestro_r44/bin/python", required = TRUE)

# organism map
org= "GRCh38" #human 
 ##mouse"GRCm38"

# run
atac <- readRDS("atac_raw_dgCMatrix.rds")
gas  <- ATACCalculateGenescore(atac, organism = org)
saveRDS(gas, "maestro.rds")
```
