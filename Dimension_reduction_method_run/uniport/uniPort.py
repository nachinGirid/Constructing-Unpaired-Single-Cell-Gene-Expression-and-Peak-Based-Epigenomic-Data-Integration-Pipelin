#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scanpy as sc
import uniport as up

### For example pbmc10k data
ds       = "pbmc"
gas_name = "signac"

# input files (user should edit paths)
rna_path = f"{ds}_rna.h5ad"
gas_path = f"{ds}_{gas_name}.h5ad"

rna = up.load_file(rna_path)
gas = up.load_file(gas_path)

# tag domains
gas.obs["domain_id"] = 0
gas.obs["source"] = "ATAC"
rna.obs["domain_id"] = 1
rna.obs["source"] = "RNA"

gas.obs["domain_id"] = gas.obs["domain_id"].astype("category")
rna.obs["domain_id"] = rna.obs["domain_id"].astype("category")

# align shared cells
shared = sorted(set(rna.obs_names).intersection(set(gas.obs_names)))
rna = rna[shared].copy()
gas = gas[shared].copy()

# filter
up.filter_data(gas, min_features=3, min_cells=200)
up.filter_data(rna, min_features=3, min_cells=200)

# concatenate on common genes
data_cm = gas.concatenate(rna, join="inner", batch_key="domain_id")

# preprocess concatenated data
sc.pp.normalize_total(data_cm)
sc.pp.log1p(data_cm)
sc.pp.highly_variable_genes(data_cm, n_top_genes=2000, subset=True)
up.batch_scale(data_cm)

# preprocess RNA
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes=2000, subset=True)
up.batch_scale(rna)

# preprocess GAS
sc.pp.normalize_total(gas)
sc.pp.log1p(gas)
sc.pp.highly_variable_genes(gas, n_top_genes=2000, subset=True)
up.batch_scale(gas)

# integrate
integrate_data = up.Run(
    adatas=[gas, rna],
    adata_cm=data_cm,
    lambda_s=1.0
)

# extract embeddings
EM  = integrate_data.obsm["latent"]
obs = integrate_data.obs

atac_mask = obs["source"].values == "ATAC"
rna_mask  = obs["source"].values == "RNA"

atac_em = pd.DataFrame(EM[atac_mask], index=obs.index[atac_mask])
rna_em  = pd.DataFrame(EM[rna_mask],  index=obs.index[rna_mask])

# remove uniPort-added suffixes
rna_em.index  = rna_em.index.str.replace(r"-(0|1)$", "", regex=True)
atac_em.index = atac_em.index.str.replace(r"-(0|1)$", "", regex=True)

# align row order
common = sorted(set(rna_em.index).intersection(set(atac_em.index)))
rna_em  = rna_em.loc[common]
atac_em = atac_em.loc[common]

# save
rna_em.to_csv(f"{ds}_rna_{gas_name}.csv")
atac_em.to_csv(f"{ds}_atac_{gas_name}.csv")
