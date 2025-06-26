import scanpy as sc

# Load filtered AnnData
adata = sc.read_h5ad('./output/cell_annotation/adata_filtered.h5ad')

# Normalize total counts per cell
sc.pp.normalize_total(adata)

# Log-transform
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata)

# Save normalized and HVG annotated AnnData
adata.write('./output/cell_annotation/adata_normalized.h5ad')
