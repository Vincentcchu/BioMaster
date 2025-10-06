#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_2
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_2/annotate_cells.py
import scanpy as sc

# Read the AnnData object with computed clusters and UMAP coordinates
adata = sc.read_h5ad("./output/cell_annotation_2/clusters_dataset.h5ad")

# Perform differential expression analysis between clusters using t-test
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test', n_genes=25)

# Annotate cells based on cluster identity (dummy criteria: cluster '0' as malignant, others as non_malignant)
adata.obs['cell_type'] = adata.obs['leiden'].apply(lambda x: "malignant" if x == "0" else "non_malignant")

# Save the annotated AnnData object
adata.write("./output/cell_annotation_2/annotated_cells.h5ad")
EOF
python ./output/cell_annotation_2/annotate_cells.py
