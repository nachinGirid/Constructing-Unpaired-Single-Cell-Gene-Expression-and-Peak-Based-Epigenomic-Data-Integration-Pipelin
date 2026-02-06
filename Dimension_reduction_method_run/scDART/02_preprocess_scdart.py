#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

### For example pbmc10k data
ds = "pbmc"
n_hvg = 1000

# input 10x folder
mtx_dir = "filtered_feature_bc_matrix"

# input region2gene
r2g_path = f"{ds}_region2gene.txt"

# output prefix
out_prefix = ds

# ---- read 10x ----
ad_all = sc.read_10x_mtx(mtx_dir, gex_only=False)

rna_vars  = ad_all.var.index[ad_all.var["feature_types"] == "Gene Expression"]
atac_vars = ad_all.var.index[ad_all.var["feature_types"] == "Peaks"]

adata_rna  = ad_all[:, rna_vars].copy()
adata_atac = ad_all[:, atac_vars].copy()

# ---- RNA HVGs ----
sc.pp.normalize_total(adata_rna, target_sum=1e4)
sc.pp.log1p(adata_rna)
sc.pp.highly_variable_genes(
    adata_rna, n_top_genes=n_hvg, flavor="seurat_v3"
)

hvg = adata_rna.var_names[adata_rna.var["highly_variable"]]

# ---- region2gene ----
r2g = pd.read_csv(
    r2g_path, sep=r"\s+", header=None,
    names=["gene", "chr", "start", "end"], engine="python"
)

r2g["peak"] = (
    r2g["chr"].astype(str) + ":" +
    r2g["start"].astype(str) + "-" +
    r2g["end"].astype(str)
)

r2g = r2g.loc[:, ["peak", "gene"]].drop_duplicates()
r2g = r2g[r2g["gene"].isin(hvg)]

# ---- subset matrices ----
peaks_keep = pd.Index(r2g["peak"]).intersection(adata_atac.var_names)
genes_keep = pd.Index(r2g["gene"]).intersection(adata_rna.var_names)

adata_atac = adata_atac[:, peaks_keep].copy()
adata_rna  = adata_rna[:,  genes_keep].copy()

def to_df(adata):
    X = adata.X
    if sparse.issparse(X):
        return pd.DataFrame.sparse.from_spmatrix(
            X, index=adata.obs_names, columns=adata.var_names
        )
    return pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)

to_df(adata_rna).to_csv(f"{out_prefix}_rna_counts.csv.gz", compression="gzip")
to_df(adata_atac).to_csv(f"{out_prefix}_atac_counts.csv.gz", compression="gzip")
r2g.to_csv(f"{out_prefix}_region2gene_hvg.csv.gz", index=False, compression="gzip")
