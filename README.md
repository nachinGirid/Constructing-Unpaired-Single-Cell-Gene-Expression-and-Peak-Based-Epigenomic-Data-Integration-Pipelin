# Benchmark_scRNA_scATAC_integration_methods
Methods are benchmarked in 3 general big steps
- Process before dimension redction
- Dimension reduction
- Labelling
## Evaluation of conversion of chromatin accessibility to gene expression: Gene activity score
four method:
1. signac
2. LIGER
3. Cicero
4. MAESTRO
### evaluation measurements
1. Correlation to the RNA measurament
2. PCA the GAS, SI of cell types
3. Anchor finding in seurat
 
 ## Evaluation of common space 
 1. celltype SI
 2. celltype cell to self distance
 3. celltype findself in 100NN
 
 ## Evaluate labelling 
 ### Is seuart cca anchor is better?
     1.cca 9543*30 ---> LIGER jointly label compare with cca tranfer label
    NMI percell==cell RNA==ATAC ATAC==RNA
     2. LIGER dimension reduction 30 ----> seurat rna clustering & transfer label
 ### IS scjoint joint clustering in Neural Network is better?
      scjoint dimension set to 30, give it to other 
 ### Is SCOT transport beter? haven't done yet 
 
 ## Methods Uses additional data that is publicly availiable
 ### GLUE
 1. celltype SI
 2. celltype cell to self distance
 3. celltype findself in 100NN
 4. with same data is this better?
 
 ### coupledNMF
 1. celltype SI
 2. celltype cell to self distance
 3. celltype findself in 100NN
 4. with same data is this better?
 
