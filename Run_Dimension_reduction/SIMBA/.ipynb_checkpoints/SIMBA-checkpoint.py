#!/usr/bin/env python3

import time
import warnings
import numpy as np
import anndata as ad
import scipy.sparse as sp
import simba as si

warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    message=r"Importing read_.* from `anndata` is deprecated"
)

### For example pbmc10k data
ds       = "pbmc"
gas_name = "signac"

start_time = time.time()

# input files (user should edit paths)
adata_CG      = ad.read_h5ad(f"{ds}_rna.h5ad")          # RNA
adata_CP      = ad.read_h5ad(f"{ds}_atac.h5ad")         # ATAC peaks
adata_CG_atac = ad.read_h5ad(f"{ds}_{gas_name}.h5ad")   # GAS

# --------- align cells across ATAC, RNA and GAS ----------
adata_CP.obs.index      = adata_CP.obs.index.astype(str)
adata_CG.obs.index      = adata_CG.obs.index.astype(str)
adata_CG_atac.obs.index = adata_CG_atac.obs.index.astype(str)

common_cells = sorted(
    set(adata_CP.obs.index)
    & set(adata_CG.obs.index)
    & set(adata_CG_atac.obs.index)
)

adata_CP      = adata_CP[common_cells].copy()
adata_CG      = adata_CG[common_cells].copy()
adata_CG_atac = adata_CG_atac[common_cells].copy()
adata_CG_atac.var_names_make_unique()

# --------- add modality suffixes ----------
adata_CP.obs.index      = adata_CP.obs.index + "_atac"
adata_CG.obs.index      = adata_CG.obs.index + "_rna"
adata_CG_atac.obs.index = adata_CP.obs.index.copy()

# --------- peak metadata ----------
adata_CP.var["gene_ids"]      = adata_CP.var.index
adata_CP.var["feature_types"] = "Peaks"
adata_CP.var["genome"]        = "hg38"

spl = np.vstack(adata_CP.var.index.astype(str).str.split("[:-]"))
adata_CP.var["chr"]   = spl[:, 0]
adata_CP.var["start"] = spl[:, 1]
adata_CP.var["end"]   = spl[:, 2]

# --------- GAS matrix sanity ----------
if "logcounts" in adata_CG_atac.layers:
    if adata_CG_atac.X.sum() == 0 and adata_CG_atac.layers["logcounts"].sum() > 0:
        adata_CG_atac.X = adata_CG_atac.layers["logcounts"].copy()

if "symbol" not in adata_CG_atac.var.columns:
    adata_CG_atac.var["symbol"] = adata_CG_atac.var.index

# ---------------- ATAC preprocessing ----------------
si.pp.filter_peaks(adata_CP, min_n_cells=3)
si.pp.cal_qc_atac(adata_CP)
si.pp.pca(adata_CP, n_components=50)
si.pp.select_pcs_features(adata_CP)

# ---------------- RNA preprocessing -----------------
si.pp.filter_genes(adata_CG, min_n_cells=3)
si.pp.cal_qc_rna(adata_CG)
si.pp.normalize(adata_CG, method="lib_size")
si.pp.log_transform(adata_CG)
si.pp.select_variable_genes(adata_CG, n_top_genes=4000)

# ---------------- GAS preprocessing -----------------
si.pp.filter_genes(adata_CG_atac, min_n_cells=3)
si.pp.cal_qc_rna(adata_CG_atac)
si.pp.normalize(adata_CG_atac, method="lib_size")
si.pp.log_transform(adata_CG_atac)

# --------- remove cells with zero shared-HVG signal ----------
def nnz_rows(X):
    if sp.issparse(X):
        return np.asarray((X > 0).sum(axis=1)).ravel()
    return (X > 0).sum(axis=1)

hvg_names = adata_CG.var_names[adata_CG.var["highly_variable"]]
shared = hvg_names.intersection(adata_CG_atac.var_names)

Xr = adata_CG[:, shared].X
Xa = adata_CG_atac[:, shared].X

good = (nnz_rows(Xr) > 0) & (nnz_rows(Xa) > 0)

adata_CG      = adata_CG[good].copy()
adata_CG_atac = adata_CG_atac[good].copy()

# --------- infer cross-modality edges ----------
adata_CrnaCatac = si.tl.infer_edges(
    adata_CG,
    adata_CG_atac,
    n_components=15,
    k=15
)
si.tl.trim_edges(adata_CrnaCatac, cutoff=0.5)

# discretize RNA
si.tl.discretize(adata_CG, n_bins=5)

# --------- build graph ----------
si.tl.gen_graph(
    list_CP=[adata_CP],
    list_CG=[adata_CG],
    list_CC=[adata_CrnaCatac],
    use_highly_variable=True,
    use_top_pcs=True,
    dirname="graph0"
)

# --------- train ----------
si.tl.pbg_train(auto_wd=True, save_wd=True, output="model")

# --------- export embeddings ----------
dict_adata = si.read_embedding()

dict_adata["C"].to_df().to_csv(
    f"{ds}_atac_{gas_name}.csv",
    sep="\t"
)

dict_adata["C2"].to_df().to_csv(
    f"{ds}_rna_{gas_name}.csv",
    sep="\t"
)

print("Elapsed time:", int(time.time() - start_time), "seconds")
