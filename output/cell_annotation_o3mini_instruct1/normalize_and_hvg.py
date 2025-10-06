import scanpy as sc

# Read the quality controlled AnnData object
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct1/adata_qc.h5ad")

# Normalize the total counts per cell to 1e4
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Save the normalized AnnData object with HVGs marked
adata.write("./output/cell_annotation_o3mini_instruct1/adata_norm.h5ad")
