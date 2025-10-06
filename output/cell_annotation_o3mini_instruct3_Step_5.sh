#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct3
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct3/cell_marker_identification.py
import scanpy as sc
import pandas as pd

# Read the AnnData object with clustering results
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct3/clustering_results.h5ad")

# Perform differential expression analysis using t-test across clusters defined by "leiden"
sc.tl.rank_genes_groups(adata, groupby="leiden", method="t-test")

# For each cluster, extract DEGs and filter for marker genes with adjusted p-value < 0.05 and log2 fold-change > 1
marker_dict = {}
# Get unique cluster labels
clusters = adata.obs["leiden"].unique().tolist()

for cluster in clusters:
    # Extract the differential expression results as a DataFrame for the given cluster
    df = sc.get.rank_genes_groups_df(adata, group=cluster)
    # Rename the columns as specified
    df = df.rename(columns={"pvals_adj": "adj_pvals", "logfoldchanges": "logFC", "pvals": "pvals"})
    # Filter for high-confidence marker genes
    df_filtered = df[(df["adj_pvals"] < 0.05) & (df["logFC"] > 1)]
    marker_dict[cluster] = df_filtered

# Save the marker gene results in the AnnData object uns attribute
adata.uns["marker_genes"] = marker_dict

# Write the updated AnnData object
adata.write_h5ad("./output/cell_annotation_o3mini_instruct3/marker_gene_results.h5ad")
EOF
python ./output/cell_annotation_o3mini_instruct3/cell_marker_identification.py
