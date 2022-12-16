# Benchmark_scRNA_scATAC_integration_methods
Methods are benchmarked in 3 general big steps
- Gene activity estimation 
- Dimension reduction
- Clustering

Three paired multiomic data is used to evaluate the performance of the Methods at different stages. At each stage, conduct the available approaches, evaluate the results using certain evaluation matrices for that stage.
## Evaluation of Gene activity score
### There are 5 options to calculate gene activity:
1. signac
2. liger
3. Cicero
4. MAESTRO
5. Unpair_reg
### evaluation measurements
1. Correlation to the RNA measurament
2. PCA the GAS, SI of cell types
3. RMSE
 
 ## Evaluation of common space 
 ### There are
 ### evaluation measurements
 1. celltype SI
 2. Distribution of distance to self
 
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
 
