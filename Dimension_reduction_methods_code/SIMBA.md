## Dimension reduction code of SIMBA  
import packages
```python
import os
import anndata as ad
import simba as si
import numpy as np
si.__version__
workdir = 'multiome_10xpmbc10k_integration'
si.settings.set_workdir(workdir)
```
Unlike the tutorial, we just read in the RNA and ATAC
```python
# subject to change for different data
atac_path="/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/seurat_object_label/atac.h5ad"
rna_path="/data2/duren_lab/naqing/pipeline_building/Data_sets/pbmc/seurat_object_label/rna.h5ad"
#Read in data
adata_CP=ad.read_h5ad(atac_path)
adata_CG=ad.read_h5ad(rna_path)
```
Operation 
```python
#peak information
adata_CP.var['gene_ids']=adata_CP.var.index
adata_CP.var['feature_types']='Peaks'
adata_CP.var['genome']='GRCh38'
n=np.array(adata_CP.var.index.str.split('[:-]'))
adata_CP.var['chr']=np.vstack(n)[:,0]
adata_CP.var['start']=np.vstack(n)[:,1]
adata_CP.var['end']=np.vstack(n)[:,2]

# manually split cells into RNAseq cells and ATACseq cells
adata_CP.obs.index = adata_CP.obs.index + '_atac'
adata_CG.obs.index = adata_CG.obs.index + '_rna'
```
Filter
```python
#ATAC
#filter peaks
si.pp.filter_peaks(adata_CP,min_n_cells=3)
#filter cells
si.pp.cal_qc_atac(adata_CP)
#si.pp.filter_cells_atac(adata_CP,min_n_peaks=100)
si.pp.pca(adata_CP, n_components=50)
si.pp.select_pcs_features(adata_CP) # didn't quite understand what is this process doing, why it is nessecary

#RNA
si.pp.filter_genes(adata_CG,min_n_cells=3)
si.pp.cal_qc_rna(adata_CG)
si.pp.normalize(adata_CG,method='lib_size')
si.pp.log_transform(adata_CG)
si.pp.select_variable_genes(adata_CG, n_top_genes=4000)
```
Discretize
```python
si.tl.discretize(adata_CG,n_bins=5) # why discretize?
```
Infer gene activity score
```python
#infer gene activity score (genome information alert!!!)
adata_CG_atac = si.tl.gene_scores(adata_CP,genome='hg38',use_gene_weigt=True, use_top_pcs=True)
#filter gene activity score
si.pp.filter_genes(adata_CG_atac,min_n_cells=3)
si.pp.cal_qc_rna(adata_CG_atac)
si.pp.normalize(adata_CG_atac,method='lib_size')
si.pp.log_transform(adata_CG_atac)
```
Infer edge and generate graph
```python
adata_CrnaCatac = si.tl.infer_edges(adata_CG, adata_CG_atac, n_components=15, k=15) # find edges
si.tl.trim_edges(adata_CrnaCatac, cutoff=0.5)# filter edges

#generate Graph
si.tl.gen_graph(list_CP=[adata_CP],
                list_CG=[adata_CG],
                list_CC=[adata_CrnaCatac],
                copy=False,
                use_highly_variable=True,
                use_top_pcs=True,
                dirname='graph0')
```                
Main step : training
```python
si.tl.pbg_train(auto_wd=True, save_wd=True, output='model')
```
Get output
```python
dict_adata=si.read_embedding() 
ATAC_em = dict_adata['C'] 
RNA_em = dict_adata['C2'] 
ATAC_em_df = ATAC_em.to_df()
RNA_em_df = RNA_em.to_df()
ATAC_em_df.to_csv("ATAC_em.txt", sep="\t", index=True)
RNA_em_df.to_csv("RNA_em.txt", sep="\t", index=True)
```

