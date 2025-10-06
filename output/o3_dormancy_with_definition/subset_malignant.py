import scanpy as sc

# Load the preprocessed annotated single-cell dataset
adata = sc.read_h5ad("./output/o3_dormancy_with_definition/preprocessed_annotated_cells.h5ad")

# Check if the 'malignant' column exists in adata.obs. If not, try using 'cell_type' for filtering
if 'malignant' in adata.obs.columns:
    malignant_cells = adata[adata.obs['malignant'] == True].copy()
elif 'cell_type' in adata.obs.columns:
    malignant_cells = adata[adata.obs['cell_type'] == 'malignant'].copy()
else:
    raise KeyError('Neither malignant nor cell_type column found in adata.obs. Please verify your dataset annotations.')

# Save the subsetted dataset
malignant_cells.write("./output/o3_dormancy_with_definition/malignant_cells.h5ad")
