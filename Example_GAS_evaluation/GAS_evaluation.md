## GAS evaluation  
Here, we use PBMC 10K's 3 methods' result as an example to show how the Gene activity score(GAS) calculating methods are evaluated.  
(The MAESTRO GAS file is too big to upload, so we only use the other methods to show examples.)
Suppose we already have the GASs calculated using different methods.  

**Preprocessing**  

As resulting GAS from different methods differ in gene number, we filter GASs by the common genes.  
Since we are using RNA of the same set of cells to evaluate the GAS, we also filter Exp(gene expression) by common genes  

```R
# Read in the GAS from different methods into one list
GASs={}
GASs[[1]]<-readRDS("pbmc_sig_sample.rds")
GASs[[2]]<-readRDS("pbmc_lig_sample.rds")
GASs[[3]]<-readRDS("pbmc_ci_sample.rds")
# Read in Exp
Gene_Exp<-readRDS("pbmc_Exp_sampe.rds")

# filter by common genes
source("GAS_eval.r")
lst<-set_common_gene(Gene_Exp,GASs)
# lst contain 3 GASs and gene_expression all on same set of genes
```
**Start evaluation**  

```R
Exp=lst$Gene_exp
GAS_list=lst$GASs
GAS_list[[3]]<-GAS_list[[3]]*10^5 ### as cicero GAS numbers are all small scale, we amplify to force them to have similar magnitude  with other GASs
```
### 1. Pearson correlation coefficient (PCC) between GAS and gene expression
```r
cor_tab<-matrix(0,3,dim(Exp)[1])
cor_tab_hv<-matrix(0,3,2000)  ### we also provide a selection of only highly variable genes correlation 
for(i in 1:3){
    correlation<-GAS_Exp_Corr(Exp,GAS_list[[i]])
    cor_tab[i,]<-correlation[[1]]
    cor_tab_hv[i,]<-correlation[[2]]
}
```
The resulting table's rows are different methods, columns are Pearson Correlation Coefficient for each gene.

### 2. Local neighborhood consistency (LNC) between cells in RNA and GAS profiles
```r
LNC<-matrix(0,3,dim(Exp)[2])

for(i in 1:3){
    LNC[i,]<-LNC(Exp,GAS_list[[i]],block=20)# neighborhood size id number of cell/20
}
```
The resulting table's rows are different methods, columns are local neighborhood consistency for each cell.

### 3. cell type Average silhouette width
```r
true_label<-read.table("pbmc_label.txt",header=TRUE,sep="\t")
GAS_list[[4]]<-Exp ## Here we RNA is evaluated together with GASs as a reference

ASW_GAS=matrix(0,length(unique(true_label$SS)),4)
for (i in 1:4){
ASW_GAS[,i]<-PCA_ASW(GAS_list[[i]],true_label$SS)
 }
colnames(ASW_GAS)<-c("Signac","liger","cicero","scRNA")
```
The resulting table's rows are cell type, columns are cell type average silhouette index for method.

```
