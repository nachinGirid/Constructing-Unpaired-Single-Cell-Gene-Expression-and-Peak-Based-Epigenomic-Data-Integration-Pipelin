## scJoint run in 3 step
- 1 preprocess, see preprocess.py
- 2 Edit config file, see below
- 3 Run main.py, directly run

### Config file (edited manually, not generated)

After preprocessing, we manually edited the scJoint config file
(config.py, copied from the scJoint MNIST/10x example) with:

number_of_class = number of unique RNA labels

input_size = number of features in the NPZ input

paths pointing to the generated .npz and rna_label.txt files

Example edits in config.py
```python
self.number_of_class = <n_class>
self.input_size = <input_size>

self.rna_paths = ['PATH/TO/RNA_signac.npz']
self.rna_labels = ['PATH/TO/rna_label.txt']
self.atac_paths = ['PATH/TO/ATAC_signac.npz']
self.atac_labels = []
```

All other parameters were left unchanged from the scJoint default
10x configuration. 