import scanpy as sc
import numpy as np

# Load the malignant cells dataset that contains cell cycle scores
adata = sc.read_h5ad('./output/o3_dormancy_with_definition/malignant_cells.h5ad')

# Verify that the cell_cycle_score column exists
if "cell_cycle_score" not in adata.obs.columns:
    raise KeyError("cell_cycle_score column not found in the dataset")

# Classify cells based on cell cycle score thresholds
conditions = [
    (adata.obs["cell_cycle_score"] < 0.3),
    (adata.obs["cell_cycle_score"] >= 0.3) & (adata.obs["cell_cycle_score"] < 0.7),
    (adata.obs["cell_cycle_score"] >= 0.7)
]
choices = ["Dormant", "Slow-cycling", "Proliferating"]

adata.obs["cell_cycle_class"] = np.select(conditions, choices, default="Undefined")

# Save the classified malignant cells dataset
adata.write('./output/o3_dormancy_with_definition/classified_malignant_cells.h5ad')
