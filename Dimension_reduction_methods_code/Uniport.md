## Uniport dimension reduction
```python
import uniport as up
import numpy as np
import pandas as pd
import scanpy as sc
print(up.__version__)
import os
```
Read in h5ad files, no ATAC need, use GAS directly
```python
rna = up.load_file("/data2/duren_lab/naqing/pipeline_building/Data_sets/BMMC_Neupris2021/Seurat_obj_labels/rna.h5ad")                    #######
gas= up.load_file("/data2/duren_lab/naqing/pipeline_building/Gene_Activity_Score/results/h5ad_GASs/gas.h5ad")
print(rna)
print(gas)
```
Put some useful information, we do not add cell type as in tutorial 
```python
rna.obs['domain_id'] = 1
rna.obs['domain_id'] = rna.obs['domain_id'].astype('category')
rna.obs['source'] = 'RNA'

gas.obs['domain_id'] = 0
gas.obs['domain_id'] = gas.obs['domain_id'].astype('category')
gas.obs['source'] = 'ATAC'
# check  if cells are matched
if (rna.obs_names.tolist() == gas.obs_names.tolist()):
    print("Row names are the same.")
else:
    print("Row names are different.")
```
Filter features and cells
```python
up.filter_data(gas, min_features=3, min_cells=200)
up.filter_data(rna, min_features=3, min_cells=200)
print(signac)
print(rna)
```
concatenate the data and preprocess 3 data (rna, gas, concatenated data)
```python
data_cm = gas.concatenate(rna, join='inner', batch_key='domain_id')
# sc.pp.highly_variable_genes(adata_cm, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(data_cm)
sc.pp.log1p(data_cm)
sc.pp.highly_variable_genes(data_cm, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(data_cm)
# sc.pp.scale(adata_cm)

# sc.pp.highly_variable_genes(adata_rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(rna)
# sc.pp.scale(adata_rna)

# sc.pp.highly_variable_genes(signac, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(gas)
sc.pp.log1p(gas)
sc.pp.highly_variable_genes(gas, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(gas)
# sc.pp.scale(gas)
```
Now integrate, main step!  
```python
integrate_data = up.Run(adatas=[gas,rna], adata_cm=data_cm, lambda_s=1.0)
```
save embedding  
```python
embedding=integrate_data.obsm['latent']
barcode=pd.DataFrame(integrate_data.obs_names)
EM= pd.DataFrame(embedding, index=barcode)
np.savetxt('bmmc_em.txt',np.hstack(barcode, EM), delimiter='\t', fmt='%s')
```          
