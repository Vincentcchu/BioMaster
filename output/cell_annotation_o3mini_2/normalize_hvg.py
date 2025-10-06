import scanpy as sc

# Read the filtered AnnData object
adata = sc.read_h5ad("./output/cell_annotation_2/filtered_dataset.h5ad")

# Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Save the normalized AnnData object annotated with HVGs
sc.write("./output/cell_annotation_2/normalized_dataset.h5ad", adata)
