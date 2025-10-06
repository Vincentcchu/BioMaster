#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_sonnet3.7t
conda install -y scanpy python-anndata leidenalg matplotlib
mkdir -p ./output/cell_annotation_sonnet3.7t
cat > ./output/cell_annotation_sonnet3.7t/explore_data.py << 'EOL'
import scanpy as sc
import os

# Load the h5ad file
adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad')

# Explore the dataset structure
summary = []
summary.append("=== Single-cell RNA-seq Dataset Exploration ===\n")

# Basic dimensions
summary.append(f"Dataset dimensions: {adata.shape[0]} cells × {adata.shape[1]} genes\n")

# Check available annotations
summary.append("=== Available Annotations ===\n")
if hasattr(adata, 'obs') and not adata.obs.empty:
    summary.append("Cell annotations (adata.obs):\n")
    for col in adata.obs.columns:
        unique_values = adata.obs[col].nunique()
        summary.append(f"  - {col}: {unique_values} unique values")
        # If there are categorical annotations with fewer than 20 categories, list them
        if unique_values < 20 and hasattr(adata.obs[col], 'cat') and adata.obs[col].dtype.name == 'category':
            categories = adata.obs[col].cat.categories.tolist()
            summary.append(f"    Categories: {', '.join(map(str, categories))}")
        summary.append("\n")
else:
    summary.append("No cell annotations available.\n")

# Check for gene annotations
summary.append("Gene annotations (adata.var):\n")
if hasattr(adata, 'var') and not adata.var.empty:
    for col in adata.var.columns:
        unique_values = adata.var[col].nunique()
        summary.append(f"  - {col}: {unique_values} unique values\n")
else:
    summary.append("No gene annotations available.\n")

# Check for layers
summary.append("=== Available Data Layers ===\n")
if hasattr(adata, 'layers') and adata.layers:
    summary.append(f"Layers: {list(adata.layers.keys())}\n")
else:
    summary.append("No additional layers available (only .X matrix).\n")

# Check for dimensionality reduction results
summary.append("=== Dimensionality Reductions ===\n")
if hasattr(adata, 'obsm') and adata.obsm:
    summary.append(f"Available embeddings: {list(adata.obsm.keys())}\n")
else:
    summary.append("No dimensionality reductions available.\n")

# Check for results from tools
summary.append("=== Analysis Results ===\n")
if hasattr(adata, 'uns') and adata.uns:
    summary.append(f"Available analysis results: {list(adata.uns.keys())}\n")
else:
    summary.append("No analysis results available.\n")

# Save the summary to a file
with open('./output/cell_annotation_sonnet3.7t/data_exploration_summary.txt', 'w') as f:
    f.write('\n'.join(summary))

print("Dataset exploration completed. Summary saved to './output/cell_annotation_sonnet3.7t/data_exploration_summary.txt'")
EOL
python ./output/cell_annotation_sonnet3.7t/explore_data.py
