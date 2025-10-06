#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct3
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct3/classify_cells.py
import scanpy as sc

# Load the AnnData object with clustering results
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct3/clustering_results.h5ad")

# Example classification logic:
# Here we assume that the clustering results contain a column "leiden" in adata.obs,
# and based on integration of marker gene profiles (e.g. high expression of known malignant markers),
# we classify specific clusters as malignant. In this example, we mark clusters '1' and '3' as malignant.

# Adjust the list below based on prior knowledge of malignant clusters
malignant_clusters = ['1', '3']

# Create a new column 'cell_class' based on the 'leiden' clustering results
adata.obs['cell_class'] = ["malignant" if cluster in malignant_clusters else "non-malignant" for cluster in adata.obs['leiden']]

# Save the final annotated AnnData object
adata.write("./output/cell_annotation_o3mini_instruct3/classified_cells.h5ad")
EOF
python ./output/cell_annotation_o3mini_instruct3/classify_cells.py
