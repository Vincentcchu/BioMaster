#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_dormancy_with_definition
conda install -y scanpy
cat << 'EOF' > ./output/o3_dormancy_with_definition/classify_malignant_cells.py\nimport scanpy as sc\nimport numpy as np\n\n# Load the malignant cells dataset that contains cell cycle scores\nadata = sc.read_h5ad('./output/o3_dormancy_with_definition/malignant_cells_cell_cycle_scores.h5ad')\n\n# Verify that the cell_cycle_score column exists\nif "cell_cycle_score" not in adata.obs.columns:\n    raise KeyError("cell_cycle_score column not found in the dataset")\n\n# Classify cells based on cell cycle score thresholds\nconditions = [\n    (adata.obs["cell_cycle_score"] < 0.3),\n    (adata.obs["cell_cycle_score"] >= 0.3) & (adata.obs["cell_cycle_score"] < 0.7),\n    (adata.obs["cell_cycle_score"] >= 0.7)\n]\nchoices = ["Dormant", "Slow-cycling", "Proliferating"]\n\nadata.obs["cell_cycle_class"] = np.select(conditions, choices, default="Undefined")\n\n# Save the classified malignant cells dataset\nadata.write('./output/o3_dormancy_with_definition/classified_malignant_cells.h5ad')\nEOF
python ./output/o3_dormancy_with_definition/classify_malignant_cells.py
