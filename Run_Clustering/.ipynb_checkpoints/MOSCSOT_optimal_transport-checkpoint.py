#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from moscot.problems.cross_modality import TranslationProblem
from scipy.stats import mode

# ----------------------------
# User-editable defaults (example-style)
# ----------------------------
METHOD         = "seurat"    # e.g., "seurat","liger","multimap","scjoint","uniport","bindsc","glue","couplednmf","scdart","simba"
DATASET        = "pbmc"      # example dataset
GAS            = "signac"    # example GAS
DIMENSION_TAG  = "dim15"     # for methods that use dim tags (seurat/liger)
# Base folder containing processed embeddings (edit once)
BASE_EMBED_DIR = Path("PATH/TO/processed_embeddings")   # e.g. processed_embeddings/
# Root with original dataset files needed for label/barcodes and rna/atac h5ad
DATA_ROOT      = Path("PATH/TO/processed_data")        # e.g. processed_data/

# ----------------------------
# Parameters
# ----------------------------
k_top   = 5
epsilon = 0.5e-2
alpha   = 0.7

# ----------------------------
# Helpers
# ----------------------------
def save_series(s: pd.Series, path_txt: Path):
    path_txt.parent.mkdir(parents=True, exist_ok=True)
    s.to_csv(path_txt, sep="\t", header=True, index=True)

def embedding_pair_paths(method, ds, gas, dimtag):
    if method in ["seurat", "liger"]:
        folder = BASE_EMBED_DIR / method / ds / gas / dimtag
    elif method in ["scjoint","uniport","bindsc","multimap"]:
        folder = BASE_EMBED_DIR / method / ds / gas
    elif method in ["glue","couplednmf","scdart","simba"]:
        folder = BASE_EMBED_DIR / method / ds
    else:
        raise ValueError(f"Unknown method: {method}")
    return folder, folder / "rna_em.csv", folder / "atac_em.csv"

def coerce_numeric_df(df: pd.DataFrame, name: str) -> pd.DataFrame:
    df_num = df.apply(pd.to_numeric, errors="coerce")
    bad_rows = df_num.isna().any(axis=1)
    if bad_rows.any():
        n_bad = int(bad_rows.sum())
        print(f"[clean] {name}: dropping {n_bad} rows with non-numeric values")
        df_num = df_num.loc[~bad_rows]
    bad_cols = df_num.isna().all(axis=0)
    if bad_cols.any():
        n_badc = int(bad_cols.sum())
        print(f"[clean] {name}: dropping {n_badc} all-NaN columns")
        df_num = df_num.loc[:, ~bad_cols]
    return df_num.astype(np.float32)

def transfer_labels_from_transport(T, rna_labels_ordered: np.ndarray, k: int):
    # T shape: (n_src_RNA, n_tgt_ATAC)
    if hasattr(T, "toarray"):
        T = T.toarray()
    topk_idx = np.argsort(-T, axis=0)[:k, :]            # (k, n_tgt)
    topk_labels = rna_labels_ordered[topk_idx]          # (k, n_tgt)
    m, _ = mode(topk_labels, axis=0, keepdims=False)    # (n_tgt,)
    return m.astype(object)

# ----------------------------
# Execution (single combo)
# ----------------------------
method = METHOD
ds     = DATASET
gas    = GAS
dimtag = DIMENSION_TAG

run_parts = [method, ds]
if method in ["seurat","liger"]:
    run_parts.extend([gas, dimtag])
else:
    if method in ["scjoint","uniport","bindsc","multimap"]:
        run_parts.append(gas)

run_key = "|".join(run_parts)

folder, rna_csv_path, atac_csv_path = embedding_pair_paths(method, ds, gas, dimtag)

print(f"[run] {run_key}")
print(f" emb folder: {folder}")

# Check embeddings exist
if not rna_csv_path.exists() or not atac_csv_path.exists():
    print(f"[skip missing] rna exists: {rna_csv_path.exists()}, atac exists: {atac_csv_path.exists()}")
    sys.exit(3)

# load embeddings
rna_em = pd.read_csv(rna_csv_path, index_col=0)
atac_em = pd.read_csv(atac_csv_path, index_col=0)
rna_em = coerce_numeric_df(rna_em, f"{run_key} RNA")
atac_em = coerce_numeric_df(atac_em, f"{run_key} ATAC")

# label and barcodes files (from original data_root)
label_file = DATA_ROOT / ds / "seurat_RNA_label.txt"
if not label_file.exists():
    alt = DATA_ROOT / ds / "seurat_RNA_label.csv"
    if alt.exists():
        label_file = alt
bc_file = DATA_ROOT / ds / "barcodes.tsv"

if not (label_file.exists() and bc_file.exists()):
    print(f"[skip missing label/barcodes] {run_key} -> {label_file}, {bc_file}")
    sys.exit(3)

lab_df = pd.read_csv(label_file, sep="\t", engine="python", header=None)
if lab_df.shape[1] != 1:
    lab_df = pd.read_csv(label_file, sep=r"\s+", engine="python", header=None)
if lab_df.shape[1] != 1:
    raise ValueError(f"{run_key}: label file has {lab_df.shape[1]} columns, expected 1")
labels = lab_df.iloc[:,0].astype(str).str.strip()

barcodes = pd.read_csv(bc_file, sep="\t", header=None, usecols=[0], engine="python")[0].astype(str)
rna_labels = pd.Series(labels.to_numpy(), index=barcodes, name="label")

# find common barcodes across embedding and labels
common = rna_em.index.intersection(atac_em.index).intersection(rna_labels.index)
if len(common) == 0:
    print(f"[skip alignment] {run_key} common=0")
    sys.exit(3)

rna_em_c  = rna_em.loc[common]
atac_em_c = atac_em.loc[common]

# coerce labels numeric (required for majority vote)
rna_lab_c = pd.to_numeric(rna_labels.loc[common], errors="raise").astype(np.int64)

# load full h5ad objects and subset in the same order
rna_adata  = sc.read_h5ad(DATA_ROOT / ds / "rna_fil_rowsum0.h5ad")
atac_adata = sc.read_h5ad(DATA_ROOT / ds / "atac_fil_rowsum0.h5ad")

rna_adata  = rna_adata[rna_adata.obs_names.isin(common)].copy()
atac_adata = atac_adata[atac_adata.obs_names.isin(common)].copy()

rna_adata  = rna_adata.loc[rna_em_c.index].copy()
atac_adata = atac_adata.loc[atac_em_c.index].copy()

rna_adata.obsm["joint"]  = rna_em_c.values
atac_adata.obsm["joint"] = atac_em_c.values

# ordered labels
rna_labels_ordered = rna_lab_c.reindex(rna_adata.obs_names).to_numpy(np.int64)

# ensure embeddings are finite
if not (np.isfinite(rna_adata.obsm["joint"]).all()):
    raise RuntimeError(f"{run_key}: non-finite values in RNA embedding")
if not (np.isfinite(atac_adata.obsm["joint"]).all()):
    raise RuntimeError(f"{run_key}: non-finite values in ATAC embedding")

# compute transport & transfer labels
ftp = TranslationProblem(adata_src=rna_adata, adata_tgt=atac_adata)
ftp = ftp.prepare(src_attr="joint", tgt_attr="joint", joint_attr="joint")
ftp = ftp.solve(epsilon=epsilon, alpha=alpha)
T = ftp.solutions[("src", "tgt")].transport_matrix

atac_pred = transfer_labels_from_transport(T, rna_labels_ordered, k_top)
out = pd.Series(atac_pred, index=atac_adata.obs_names, name="label")

save_series(out, folder / "atac_label_moscot.txt")
print(f"[done] {run_key}")
