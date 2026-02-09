#!/usr/bin/env python3

import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from moscot.problems.cross_modality import TranslationProblem
from scipy.stats import mode

warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", FutureWarning)

### For example pbmc10k data
method = "seurat"
ds     = "pbmc"
gas    = "signac"
dim    = "dim15"

# parameters
k_top   = 5
epsilon = 0.5e-2
alpha   = 0.7

# input folders (user should edit paths)
embed_dir = Path("PATH/TO/processed_embeddings") / method / ds / gas / dim
data_root = Path("PATH/TO/processed_data") / ds

rna_csv  = embed_dir / "rna_em.csv"
atac_csv = embed_dir / "atac_em.csv"

# --------------------
# Load embeddings
# --------------------
rna_em  = pd.read_csv(rna_csv,  index_col=0)
atac_em = pd.read_csv(atac_csv, index_col=0)

# ensure numeric
rna_em  = rna_em.apply(pd.to_numeric)
atac_em = atac_em.apply(pd.to_numeric)

# --------------------
# Load RNA labels + barcodes
# --------------------
label_file = data_root / "seurat_RNA_label.txt"
if not label_file.exists():
    label_file = data_root / "seurat_RNA_label.csv"

bc_file = data_root / "barcodes.tsv"

labels = pd.read_csv(label_file, sep=r"\s+", header=None)[0].astype(str)
barcodes = pd.read_csv(bc_file, sep="\t", header=None)[0].astype(str)

rna_labels = pd.Series(labels.to_numpy(), index=barcodes, name="label")

# --------------------
# Align shared cells
# --------------------
common = rna_em.index.intersection(atac_em.index).intersection(rna_labels.index)

rna_em  = rna_em.loc[common]
atac_em = atac_em.loc[common]
rna_lab = pd.to_numeric(rna_labels.loc[common], errors="raise").astype(np.int64)

# --------------------
# Load AnnData objects
# --------------------
rna_adata  = sc.read_h5ad(data_root / "rna_fil_rowsum0.h5ad")
atac_adata = sc.read_h5ad(data_root / "atac_fil_rowsum0.h5ad")

rna_adata  = rna_adata[rna_adata.obs_names.isin(common)].copy()
atac_adata = atac_adata[atac_adata.obs_names.isin(common)].copy()

rna_adata  = rna_adata.loc[rna_em.index].copy()
atac_adata = atac_adata.loc[atac_em.index].copy()

rna_adata.obsm["joint"]  = rna_em.values.astype(np.float32)
atac_adata.obsm["joint"] = atac_em.values.astype(np.float32)

rna_labels_ordered = rna_lab.reindex(rna_adata.obs_names).to_numpy()

# --------------------
# MOSCOT label transfer
# --------------------
tp = TranslationProblem(
    adata_src=rna_adata,
    adata_tgt=atac_adata
)

tp = tp.prepare(
    src_attr="joint",
    tgt_attr="joint",
    joint_attr="joint"
)

tp = tp.solve(
    epsilon=epsilon,
    alpha=alpha
)

T = tp.solutions[("src", "tgt")].transport_matrix
if hasattr(T, "toarray"):
    T = T.toarray()

# top-k majority vote
topk_idx = np.argsort(-T, axis=0)[:k_top, :]
topk_lab = rna_labels_ordered[topk_idx]
pred, _  = mode(topk_lab, axis=0, keepdims=False)

# --------------------
# Save output
# --------------------
out = pd.Series(
    pred.astype(str),
    index=atac_adata.obs_names,
    name="label"
)

out.to_csv(embed_dir / "atac_label_moscot.txt", sep="\t")

print("Done: MOSCOT label transfer")
