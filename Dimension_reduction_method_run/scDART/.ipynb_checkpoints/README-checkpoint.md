## scDART (example: pbmc)

This folder contains a minimal, example-based implementation of scDART,
split into three sequential steps.

### Step 1: region-to-gene mapping
Generate a long-format region–gene table using bedtools:

```bash
bash 01_make_region2gene.sh
```
Output:
- pbmc_region2gene.txt

### Step 2: preprocessing

Select RNA HVGs, filter ATAC peaks and region–gene pairs, and export
count matrices:

```bash
python 02_preprocess_scdart.py
```
Outputs:

- pbmc_rna_counts.csv.gz
- pbmc_atac_counts.csv.gz
- pbmc_region2gene_hvg.csv.gz

### Step 3: scDART integration

Run scDART using parameters yield lowest MMD loss as used in the benchmark:
```bash
python 03_run_scdart.py
```

Outputs:
- pbmc_z_rna.csv.gz
- pbmc_z_atac.csv.gz
- pbmc_scdart_meta.json

All scripts use pbmc as an example dataset and are intended as
reference implementations rather than turnkey pipelines.