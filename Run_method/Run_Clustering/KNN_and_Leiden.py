#!/usr/bin/env python3

import warnings
from pathlib import Path
import pandas as pd
import scanpy as sc
from sklearn.neighbors import KNeighborsClassifier

warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", FutureWarning)

### For example pbmc10k data
method = "seurat"
ds     = "pbmc"
gas    = "signac"
dim    = "dim15"

# parameters
k_neighbors        = 5
leiden_res         = 1.0
n_neighbors_graph  = 20
metric             = "cosine"

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
rna_lab = rna_labels.loc[common]

# --------------------
# KNN label transfer
# --------------------
knn = KNeighborsClassifier(n_neighbors=k_neighbors)
knn.fit(rna_em.values, rna_lab.values)
atac_pred = knn.predict(atac_em.values)

pd.Series(
    atac_pred,
    index=atac_em.index,
    name="label"
).to_csv(embed_dir / "atac_label_knn.txt", sep="\t")

# --------------------
# Leiden clustering (joint)
# --------------------
rna_u  = rna_em.copy();  rna_u.index  = rna_u.index.astype(str)  + "_RNA"
atac_u = atac_em.copy(); atac_u.index = atac_u.index.astype(str) + "_ATAC"

combined = pd.concat([rna_u, atac_u], axis=0)

ad = sc.AnnData(combined.values)
ad.obs_names = combined.index

sc.pp.neighbors(ad, n_neighbors=n_neighbors_graph, metric=metric)
sc.tl.leiden(ad, resolution=leiden_res)

clusters = ad.obs["leiden"].copy()
clusters.name = "cluster"

rna_cl  = clusters[clusters.index.str.endswith("_RNA")]
atac_cl = clusters[clusters.index.str.endswith("_ATAC")]

rna_cl.index  = rna_cl.index.str.replace("_RNA$",  "", regex=True)
atac_cl.index = atac_cl.index.str.replace("_ATAC$", "", regex=True)

rna_cl.to_csv(embed_dir / "rna_label_leiden.txt",  sep="\t")
atac_cl.to_csv(embed_dir / "atac_label_leiden.txt", sep="\t")

print("Done: KNN label transfer and Leiden clustering")
