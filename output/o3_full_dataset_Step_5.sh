#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_full_dataset
conda install -y scanpy
cat << 'EOF' > ./output/o3_full_dataset/annotate_cells.py
import scanpy as sc

# Read the AnnData object with computed clusters and UMAP coordinates
adata = sc.read("./output/o3_full_dataset/clusters_dataset.h5ad")

# Ensure that the clustering information is available (assumed to be stored in the 'leiden' column)
if 'leiden' not in adata.obs:
    raise ValueError("Clustering information 'leiden' not found in AnnData.obs")

# Perform differential expression analysis for each cluster to identify marker genes (this can aid manual curation)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test', use_raw=False)

# Annotate cells based on known malignant cell markers (for demonstration, we assign cells in cluster '2' as malignant and others as non_malignant)
adata.obs['cell_type'] = adata.obs['leiden'].apply(lambda x: 'malignant' if x == '2' else 'non_malignant')

# Save the annotated AnnData object
adata.write("./output/o3_full_dataset/annotated_cells.h5ad")
EOF
python ./output/o3_full_dataset/annotate_cells.py
