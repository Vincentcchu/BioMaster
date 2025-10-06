import scanpy as sc

# Read the preprocessed AnnData object with quality controlled and normalized data
adata = sc.read_h5ad('./output/o3_dormancy_with_malignancy2/preprocessed_data.h5ad')

# Filter the cells to retain only those labeled as 'malignant'
# Assuming that the malignancy annotations are stored in the adata.obs column named 'malignancy'
malignant_adata = adata[adata.obs['malignancy'] == 'malignant'].copy()

# Write the filtered AnnData object to the output file
malignant_adata.write('./output/o3_dormancy_with_malignancy2/malignant_cells.h5ad')
