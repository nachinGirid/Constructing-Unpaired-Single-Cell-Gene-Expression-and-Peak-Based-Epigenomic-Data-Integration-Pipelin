#!/usr/bin/env python3

import time
import pandas as pd
import scanpy as sc
import anndata as ad
import scglue
from itertools import chain

### For example pbmc10k data
ds = "pbmc"

# input files (user should edit paths)
rna = ad.read_h5ad(f"{ds}_rna.h5ad")
atac = ad.read_h5ad(f"{ds}_atac.h5ad")

# gene annotation (example: mouse mm10)
gtf = "PATH/TO/GTF"

# --------------------
# RNA preprocessing
# --------------------
rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100)

# --------------------
# ATAC preprocessing
# --------------------
scglue.data.lsi(atac, n_components=100, n_iter=15)

# --------------------
# Gene annotation
# --------------------
scglue.data.get_gene_annotation(rna, gtf=gtf, gtf_by="gene_name")
rna = rna[:, rna.var.dropna(subset=["chrom"]).index]

# --------------------
# Peak coordinates
# --------------------
split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: int(x[1]))
atac.var["chromEnd"] = split.map(lambda x: int(x[2]))

# --------------------
# Guidance graph
# --------------------
guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
scglue.graph.check_graph(guidance, [rna, atac])

# --------------------
# Configure datasets
# --------------------
scglue.models.configure_dataset(
    rna,
    distribution="NB",
    use_highly_variable=True,
    use_layer="counts",
    use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac,
    distribution="NB",
    use_highly_variable=True,
    use_rep="X_lsi"
)

guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()


# --------------------
# Train GLUE
# --------------------
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac},
    guidance_hvf,
    fit_kws={"directory": "glue"}
)

print("Integration consistency:")
print(scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance_hvf
))

# --------------------
# Save embeddings
# --------------------
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

pd.DataFrame(
    rna.obsm["X_glue"], index=rna.obs_names
).to_csv(f"{ds}_rna_glue.csv")

pd.DataFrame(
    atac.obsm["X_glue"], index=atac.obs_names
).to_csv(f"{ds}_atac_glue.csv")
