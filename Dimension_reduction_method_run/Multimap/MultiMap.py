#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import MultiMAP
from scipy.sparse import csr_matrix

### For example pbmc10k data
ds  = "pbmc"
gas = "signac"

# input files (user should edit paths)
rna        = sc.read(f"{ds}_rna.h5ad")
atac_peaks = sc.read(f"{ds}_atac.h5ad")
atac_genes = sc.read(f"{ds}_{gas}.h5ad")

# tag barcodes by modality (treat as unpaired)
for aobj, suf in (
    (rna, "_rna"),
    (atac_peaks, "_atac"),
    (atac_genes, "_atac")
):
    aobj.obs_names = aobj.obs_names.astype(str) + suf
    aobj.obs_names_make_unique()

rna_cells  = set(rna.obs_names)
atac_cells = set(atac_genes.obs_names)

# ensure GAS object has all ATAC cells
missing = [c for c in atac_peaks.obs_names if c not in atac_genes.obs_names]
if len(missing) > 0:
    filler = ad.AnnData(
        X=csr_matrix((len(missing), atac_genes.n_vars)),
        obs=pd.DataFrame(index=missing),
        var=atac_genes.var.copy(),
    )
    atac_genes = ad.concat([atac_genes, filler], axis=0)

# match cell order
atac_genes = atac_genes[atac_peaks.obs_names, :]

# LSI on peaks, copy into GAS object
MultiMAP.TFIDF_LSI(atac_peaks)
atac_genes.obsm["X_lsi"] = atac_peaks.obsm["X_lsi"].copy()

# PCA on RNA
rna_pca = rna.copy()
sc.pp.scale(rna_pca)
sc.pp.pca(rna_pca)
rna.obsm["X_pca"] = rna_pca.obsm["X_pca"].copy()

# integrate (default parameters)
adata = MultiMAP.Integration(
    [rna, atac_genes],
    ["X_pca", "X_lsi"],
    strengths=[0.8, 0.2]
)

# modality labels
mod = np.where(
    np.isin(adata.obs_names, list(rna_cells)), "RNA",
    np.where(np.isin(adata.obs_names, list(atac_cells)), "ATAC", "NA")
)
adata.obs["modality"] = mod

# save embedding
emb = adata.obsm["X_multimap"]
emb_df = pd.DataFrame(
    emb,
    index=adata.obs_names,
    columns=[f"MM_{i+1}" for i in range(emb.shape[1])]
)
emb_df["modality"] = adata.obs["modality"].values

emb_df.to_csv(f"{ds}_{gas}_multimap.csv")

# save AnnData
adata.write(f"{ds}_{gas}_multimap.h5ad", compression="gzip")
