import scanpy as sc
import sys

# Load the AnnData object containing the raw malignant cell data
adata = sc.read_h5ad('./output/o3_dormancy/initialized_malignant_cells.h5ad')

# Load the malignant classification AnnData object
classification = sc.read_h5ad('/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/biomaster/output/cell_annotation_o3mini/malignant_classification.h5ad')

# Determine the correct column name for malignant labels
if 'malignant' in classification.obs.columns:
    col_name = 'malignant'
elif 'malignant_label' in classification.obs.columns:
    col_name = 'malignant_label'
else:
    sys.exit('Error missing malignant label in classification.obs')

# Join the classification labels to the AnnData object
adata.obs = adata.obs.join(classification.obs[[col_name]], how='left')

# Subset cells that are labeled as malignant
adata_malignant = adata[adata.obs[col_name] == True].copy()

# Write the subset AnnData object to file
sc.write(adata_malignant, './output/o3_dormancy/subset_malignant_cells.h5ad')
