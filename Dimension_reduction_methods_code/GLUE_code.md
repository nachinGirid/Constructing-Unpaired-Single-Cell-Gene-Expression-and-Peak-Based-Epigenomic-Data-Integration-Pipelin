# GLUE
## convert seurat object to h5ad
```r
library(Seurat)
library(zellkonverter)
counts <- Read10X_h5('/data2/duren_lab/naqing/data/HumanBrain/human_brain_3k_filtered_feature_bc_matrix.h5')

rna <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

atac <- CreateSeuratObject(
  counts = counts$`Peak`,
  assay = "ATAC"
)
rna<-as.SingleCellExperiment(pbmc.rna)
atac<-as.SingleCellExperiment(pbmc.atac)

writeH5AD(rna,"/data/duren_lab/naqing/Methods_Benchmark/GLUE/pbmc.rna.h5ad")
writeH5AD(atac,"/data/duren_lab/naqing/Methods_Benchmark/GLUE/pbmc.atac.h5ad")
```
## Stage 1: Data processing
### preprocess data
```python
import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)
rna = ad.read_h5ad("pbmc.rna.h5ad")
atac = ad.read_h5ad("pbmc.atac.h5ad")
# RNA
rna.X, rna.X.data
rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")
sc.pp.neighbors(rna, metric="cosine")
# ATAC
atac.X, atac.X.data
scglue.data.lsi(atac, n_components=100, n_iter=15)
```
### construct guidance graph
```python
# RNA obtain coordinate
scglue.data.get_gene_annotation(
    rna, gtf="/data/duren_lab/naqing/Methods_Benchmark/GLUE/gencode.v41.chr_patch_hapl_scaff.basic.annotation.gtf.gz",
    gtf_by="gene_name"
)
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
rna.var.dtypes
rna.var['artif_dupl']=False
# atac
split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
atac.var.head()
# build
guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
guidance
scglue.graph.check_graph(guidance, [rna, atac])
```
## Stage 2: Model training
Read in stored data if needed:
```python
rna = ad.read_h5ad("rna-pp.h5ad")
atac = ad.read_h5ad("atac-pp.h5ad")
guidance = nx.read_graphml("guidance.graphml.gz")
```
### Configure data
```python
from itertools import chain
import itertools
import networkx as nx
import pandas as pd
import seaborn as sns
# Here we can decide if use only highly variable genes, I used, because it is THE default
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)
scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)
# Since only highly variable genes are used, graph only for these genes are need 
guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()
```
### Train GLUE model
```python
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": "glue"}
)
# result saving
# glue.save("glue.dill")
```
```python
# if new session, please load the former object and libraries
glue = scglue.models.load_model("glue.dill")
glue.compile()
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance_hvf
)
dx # confidency is good 
```
### Extract the embeddings
```python
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
cell_embeddings = pd.DataFrame(rna.obsm['X_glue'])
cell_embeddings.iloc[:5, :5]
```

### evaluate he embedding
```r
source("/data2/duren_lab/naqing/pipeline_building/functions/Func_DMR_eval.r")
rna_label<-readRDS('/data2/duren_lab/naqing/Benchmark_human-brain/seurat_RNA_label.rds')
# human
setwd('/data2/duren_lab/naqing/Methods_Benchmark/GLUE/human_data/')
rna_em<-read.table('rna_embedding.txt',sep='\t')
RNA_em<-rna_em[,2:dim(rna_em)[2]]
rownames(RNA_em)<-rna_em[,1]
atac_em<-read.table('atac_embedding.txt',sep='\t')
ATAC_em<-atac_em[,2:dim(atac_em)[2]]
rownames(ATAC_em)<-atac_em[,1]
eval_res<-embedding_evaluation(RNA_em,ATAC_em,rna_label)
saveRDS(eval_res,'GLUE_DMR_eval.rds')
```
