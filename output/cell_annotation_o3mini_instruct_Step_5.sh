#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct
conda install -y -c conda-forge scanpy pandas
mkdir -p ./output/cell_annotation_o3mini_instruct
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct/marker_gene_identification.py
import scanpy as sc
import pandas as pd

# Load the normalized AnnData object
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct/normalized_hvg_dataset.h5ad")

# Load clustering information and add it to the AnnData object
adata_cluster = sc.read_h5ad("./output/cell_annotation_o3mini_instruct/clustering_results.h5ad")
adata.obs['leiden'] = adata_cluster.obs['leiden']

# Perform differential expression analysis using t-test across clusters
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')

# Extract results, rename columns, and filter based on criteria (adj_pvals < 0.05 and logFC > 1)
groups = adata.obs['leiden'].unique()
marker_list = []
for group in groups:
    df = sc.get.rank_genes_groups_df(adata, group=group)
    # Rename columns to required output names
    df = df.rename(columns={"logfoldchanges": "logFC", "pvals": "pvals", "pvals_adj": "adj_pvals"})
    # Filter for high-confidence markers
    df_filtered = df[(df["adj_pvals"] < 0.05) & (df["logFC"] > 1)]
    df_filtered['cluster'] = group
    marker_list.append(df_filtered)

# Combine results from all clusters
markers_df = pd.concat(marker_list, ignore_index=True)

# Save the marker genes table
markers_df.to_csv("./output/cell_annotation_o3mini_instruct/marker_genes_table.txt", sep="\t", index=False)
EOF
python ./output/cell_annotation_o3mini_instruct/marker_gene_identification.py
