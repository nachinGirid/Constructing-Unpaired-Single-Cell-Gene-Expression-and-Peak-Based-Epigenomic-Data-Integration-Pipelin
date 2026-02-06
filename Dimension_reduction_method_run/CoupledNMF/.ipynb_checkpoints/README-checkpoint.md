## CoupledNMF

CoupledNMF was run using the original MATLAB implementation provided by the authors.

### Running CoupledNMF

The main script is executed in MATLAB as follows:

```bash
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "code_coupleNMF; exit"
```
Input format

The input directory (Indir) provided to the script must contain a 10x Genomics–style matrix with the following files:

- matrix.mtx
- features.tsv
- barcodes.tsv

These files define a sparse count matrix in Matrix Market format, following the standard 10x Genomics convention.

**Notes**

- All parameters were left at their default values as provided in the original CoupledNMF code.
- All functions called by code_coupleNMF.m are included in the same folder as the main script.