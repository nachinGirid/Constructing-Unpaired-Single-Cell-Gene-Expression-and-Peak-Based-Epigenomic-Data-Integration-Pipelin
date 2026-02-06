#!/usr/bin/env python3

import numpy as np
import pandas as pd
import torch
import scDART
import json

### For example pbmc10k data
ds = "pbmc"

# fixed parameters (as used in benchmark)
latent_dim = 8
reg_g = 1
reg_mmd = 10
epochs = 500
ts = [30, 50, 70]
LR = 3e-4
seed = 0

np.random.seed(seed)
torch.manual_seed(seed)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# inputs
rna  = pd.read_csv(f"{ds}_rna_counts.csv.gz", index_col=0)
atac = pd.read_csv(f"{ds}_atac_counts.csv.gz", index_col=0)
reg_long = pd.read_csv(f"{ds}_region2gene_hvg.csv.gz")

# build region matrix
reg_long["value"] = 1
reg = reg_long.pivot_table(
    index="peak", columns="gene", values="value", fill_value=0
)

rna  = rna.loc[:, reg.columns]
atac = atac.loc[:, reg.index]

model = scDART.scDART(
    n_epochs=epochs,
    latent_dim=latent_dim,
    ts=ts,
    use_anchor=False,
    use_potential=True,
    k=10,
    reg_d=1,
    reg_g=reg_g,
    reg_mmd=reg_mmd,
    l_dist_type="kl",
    seed=seed,
    learning_rate=LR,
    device=device
).fit(
    rna_count=rna.values.astype(np.float32),
    atac_count=atac.values.astype(np.float32),
    reg=reg.values.astype(np.float32),
    rna_anchor=None,
    atac_anchor=None
)

z_rna, z_atac = model.transform(
    rna.values.astype(np.float32),
    atac.values.astype(np.float32)
)

cols = [f"z{i+1}" for i in range(z_rna.shape[1])]
pd.DataFrame(z_rna, index=rna.index, columns=cols)\
  .to_csv(f"{ds}_z_rna.csv.gz", compression="gzip")
pd.DataFrame(z_atac, index=atac.index, columns=cols)\
  .to_csv(f"{ds}_z_atac.csv.gz", compression="gzip")

with open(f"{ds}_scdart_meta.json", "w") as f:
    json.dump({
        "latent_dim": latent_dim,
        "reg_g": reg_g,
        "reg_mmd": reg_mmd,
        "epochs": epochs,
        "ts": ts,
        "learning_rate": LR,
        "seed": seed
    }, f, indent=2)
