#scJoint preprocessing (example: pbmc + signac)
#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import process_db   # scJoint preprocessing utilities

### For example pbmc10k data
ds  = "pbmc"
gas = "signac"

# input directory (user should edit)
basepath = "PATH/TO/SCJOINT_INPUT"

# expected inputs (prepared earlier)
rna_h5  = f"{basepath}/{ds}/{gas}/RNA_{gas}.h5"
atac_h5 = f"{basepath}/{ds}/{gas}/ATAC_{gas}.h5"
rna_csv = f"{basepath}/{ds}/{gas}/rna_label.csv"

# sanity check
for p in [rna_h5, atac_h5, rna_csv]:
    if not os.path.exists(p):
        raise FileNotFoundError(p)

# scJoint expects lists of paths
rna_h5_files     = [rna_h5]
atac_h5_files    = [atac_h5]
rna_label_files  = [rna_csv]
atac_label_files = []

# --------------------
# Step 1: data parsing
# --------------------
# produces:
#   RNA_<gas>.npz
#   ATAC_<gas>.npz
process_db.data_parsing(rna_h5_files, atac_h5_files)

# --------------------
# Step 2: label parsing
# --------------------
# produces:
#   rna_label.txt (numeric labels)
process_db.label_parsing(rna_label_files, atac_label_files)

# --------------------
# Step 3: inspect outputs
# --------------------
label_txt = f"{basepath}/{ds}/{gas}/rna_label.txt"
labels = pd.read_csv(label_txt, header=None, sep="\t")
n_class = labels[0].nunique()

rna_npz = f"{basepath}/{ds}/{gas}/RNA_{gas}.npz"
with np.load(rna_npz) as npz:
    input_size = int(npz["shape"][1])

print("scJoint preprocessing finished")
print(f"number_of_class = {n_class}")
print(f"input_size     = {input_size}")
