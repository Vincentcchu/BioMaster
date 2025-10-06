#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct2
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct2/marker_gene_identification.py
import scanpy as sc
import pandas as pd

# Load the normalized AnnData object
adata = sc.read_h5ad('./output/cell_annotation_o3mini_instruct2/normalized_dataset.h5ad')

# If the leiden clustering labels are missing then compute neighbors and perform leiden clustering
if 'leiden' not in adata.obs.columns:
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

# Perform differential expression analysis using t-test
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')

# Collect marker genes for each cluster
marker_list = []
for group in adata.obs['leiden'].unique():
    df = sc.get.rank_genes_groups_df(adata, group=group)
    df.rename(columns={"logfoldchanges": "logFC", "pvals_adj": "adj_pvals"}, inplace=True)
    df_filtered = df[(df['adj_pvals'] < 0.05) & (df['logFC'] > 1)]
    df_filtered['cluster'] = group
    marker_list.append(df_filtered)

# Concatenate results and save to file
marker_genes = pd.concat(marker_list, ignore_index=True)
marker_genes.to_csv('./output/cell_annotation_o3mini_instruct2/marker_genes.txt', sep='\t', index=False)

print("Marker gene identification completed. Results saved to ./output/cell_annotation_o3mini_instruct2/marker_genes.txt")
EOF
python ./output/cell_annotation_o3mini_instruct2/marker_gene_identification.py
