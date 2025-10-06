#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_sonnet3.7t
conda install -y scanpy pandas numpy matplotlib
cat << 'EOF' > ./output/cell_annotation_sonnet3.7t/identify_marker_genes.py
# Import necessary libraries
import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Read the dimensionality reduction and clustering results
print("Reading data...")
adata = sc.read_h5ad('./output/cell_annotation_sonnet3.7t/dimensionality_reduction_clustering.h5ad')

# Check if clustering has been performed
if 'leiden' not in adata.obs.columns:
    print("Leiden clustering not found. Checking for other clustering results...")
    cluster_keys = [key for key in adata.obs.columns if 'cluster' in key.lower() or 'louvain' in key.lower()]
    if cluster_keys:
        cluster_key = cluster_keys[0]
        print(f"Using {cluster_key} for cluster identification.")
    else:
        raise ValueError("No clustering information found in the data. Please perform clustering first.")
else:
    cluster_key = 'leiden'
    print("Using Leiden clustering for marker gene identification.")

# Identify marker genes for each cluster
print(f"Identifying marker genes for each {cluster_key} cluster...")
try:
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='t-test')
except Exception as e:
    print(f"Error in rank_genes_groups: {e}")
    print("Trying with different parameters...")
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon')

# Extract the results and convert to a DataFrame
marker_genes_dfs = []
clusters = adata.obs[cluster_key].unique()
print(f"Found {len(clusters)} clusters.")

for cluster in clusters:
    try:
        cluster_str = str(cluster)
        gene_names = adata.uns['rank_genes_groups']['names'][cluster_str]
        logfcs = adata.uns['rank_genes_groups']['logfoldchanges'][cluster_str]
        pvals = adata.uns['rank_genes_groups']['pvals'][cluster_str]
        pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][cluster_str]
        
        genes = pd.DataFrame({
            'gene': gene_names,
            'logFC': logfcs,
            'pvals': pvals,
            'adj_pvals': pvals_adj,
            'cluster': cluster
        })
        
        # Filter for significant genes
        significant_genes = genes[(genes['adj_pvals'] < 0.05) & (genes['logFC'] > 1)]
        print(f"Cluster {cluster}: {len(significant_genes)} significant marker genes.")
        marker_genes_dfs.append(significant_genes)
    except Exception as e:
        print(f"Error extracting marker genes for cluster {cluster}: {e}")
        continue

# Combine all marker genes into a single DataFrame
if marker_genes_dfs:
    all_markers = pd.concat(marker_genes_dfs)
    print(f"Found a total of {len(all_markers)} marker genes across all clusters.")
    
    # Save marker genes to CSV
    all_markers.to_csv('./output/cell_annotation_sonnet3.7t/marker_genes.csv', index=False)
    print("Marker genes saved to './output/cell_annotation_sonnet3.7t/marker_genes.csv'")
else:
    print("No marker genes found. Check the clustering results and parameters.")
    all_markers = pd.DataFrame(columns=['gene', 'logFC', 'pvals', 'adj_pvals', 'cluster'])
    all_markers.to_csv('./output/cell_annotation_sonnet3.7t/marker_genes.csv', index=False)

# Define a list of genes commonly associated with malignancy
malignancy_genes = [
    'MYC', 'EGFR', 'KRAS', 'TP53', 'BRCA1', 'BRCA2', 'PTEN', 'RB1', 
    'CDKN2A', 'PIK3CA', 'CCND1', 'CDK4', 'MDM2', 'ERBB2', 'VEGFA',
    'MMP2', 'MMP9', 'BIRC5', 'BCL2', 'BAX', 'CASP3', 'CASP9',
    'CDH1', 'CTNNB1', 'SNAI1', 'SNAI2', 'ZEB1', 'ZEB2', 'TWIST1',
    'IL6', 'IL8', 'TNF', 'TGFB1', 'IFNG', 'STAT3', 'NFKB1',
    'SOX2', 'POU5F1', 'NANOG', 'CD44', 'PROM1', 'ALDH1A1'
]

# Identify which marker genes are associated with malignancy
if not all_markers.empty:
    malignancy_markers = all_markers[all_markers['gene'].isin(malignancy_genes)]
    print(f"Found {len(malignancy_markers)} marker genes associated with malignancy.")
    
    # Save malignancy gene signature to CSV
    malignancy_markers.to_csv('./output/cell_annotation_sonnet3.7t/malignancy_gene_signature.csv', index=False)
    print("Malignancy gene signature saved to './output/cell_annotation_sonnet3.7t/malignancy_gene_signature.csv'")
else:
    print("No marker genes found, so no malignancy gene signature could be generated.")
    malignancy_markers = pd.DataFrame(columns=['gene', 'logFC', 'pvals', 'adj_pvals', 'cluster'])
    malignancy_markers.to_csv('./output/cell_annotation_sonnet3.7t/malignancy_gene_signature.csv', index=False)

# Visualize top marker genes per cluster
if not all_markers.empty:
    try:
        top_markers = all_markers.groupby('cluster').apply(lambda x: x.nlargest(10, 'logFC')).reset_index(drop=True)
        sc.pl.heatmap(adata, top_markers['gene'].unique(), groupby=cluster_key, 
                      dendrogram=True, standard_scale='var', swap_axes=True,
                      save='_top_marker_genes.pdf')
    except Exception as e:
        print(f"Error generating heatmap: {e}")
EOF
python ./output/cell_annotation_sonnet3.7t/identify_marker_genes.py
