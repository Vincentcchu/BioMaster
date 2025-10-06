#!/usr/bin/env python
# Script to perform dimensionality reduction and clustering on single-cell data

import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import os
import numpy as np

# Ensure output directory exists
os.makedirs('./output/cell_annotation_sonnet3.7t/', exist_ok=True)

# Configure matplotlib manually
plt.rcParams['figure.figsize'] = (10, 10)
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 150

# Set basic scanpy settings - avoid the problematic set_figure_params
sc.settings.verbosity = 3
sc.settings.figdir = './output/cell_annotation_sonnet3.7t/'

# Read preprocessed data
print("Reading preprocessed data...")
adata = sc.read_h5ad('./output/cell_annotation_sonnet3.7t/preprocessed_data.h5ad')
print(f"Data shape: {adata.shape}")

# Perform PCA for dimensionality reduction
print("Performing PCA...")
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# Compute neighborhood graph
print("Computing neighborhood graph...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# Perform UMAP for visualization
print("Performing UMAP...")
sc.tl.umap(adata)

# Apply clustering to identify cell populations
print("Performing Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5)

# Save the anndata object with dimensionality reduction and clustering results
print("Saving processed data...")
adata.write('./output/cell_annotation_sonnet3.7t/dimensionality_reduction_clustering.h5ad')

# Create visualizations using the most robust approach
print("Creating visualizations...")

# Create UMAP visualization with leiden clusters
fig, ax = plt.subplots(figsize=(10, 10))
sc.pl.umap(adata, color=['leiden'], title='Cell Clusters (Leiden)', ax=ax, show=False)
fig.savefig('./output/cell_annotation_sonnet3.7t/clustering_visualization.pdf', bbox_inches='tight')
plt.close(fig)

print("Analysis complete!")
