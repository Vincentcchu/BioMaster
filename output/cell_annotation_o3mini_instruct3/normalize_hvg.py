import scanpy as sc

# Read the QC-filtered AnnData object
adata = sc.read_h5ad('./output/cell_annotation_o3mini_instruct3/qc_filtered_data.h5ad')

# Normalize the data by total counts per cell
sc.pp.normalize_total(adata)

# Log-transform the data
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata)

# Save the normalized data with HVG information
adata.write('./output/cell_annotation_o3mini_instruct3/normalized_hvg_data.h5ad')
