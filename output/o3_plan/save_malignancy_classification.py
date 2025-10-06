import scanpy as sc

# Load the AnnData object with malignancy scores and binary classification labels
adata = sc.read_h5ad('./output/o3_plan/scored_classified_dataset.h5ad')

# Save the final AnnData object with malignancy annotations
adata.write('./output/o3_plan/malignancy_classification.h5ad')
