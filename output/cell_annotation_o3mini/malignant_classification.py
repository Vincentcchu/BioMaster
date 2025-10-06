#!/usr/bin/env python

import scanpy as sc
import celltypist
import pandas as pd

# Load the AnnData object with cluster labels (assumes dimensionality reduction and clustering results are stored here)
adata = sc.read_h5ad('./output/cell_annotation/clusters_annotated.h5ad')

# Load the differential expression marker gene list for clusters
markers = pd.read_csv('./output/cell_annotation/cluster_markers.csv')

# Define known malignant markers (update this list with current literature markers as needed)
malignant_markers = ['GeneA', 'GeneB']

# Example strategy: assign malignant classification based on cluster label
# Here, we assume that cluster '0' is malignant based on marker analysis
adata.obs['cell_class'] = adata.obs['leiden'].apply(lambda x: 'malignant' if x == '0' else 'non-malignant')

# Optionally, one could integrate celltypist predictions or use expression thresholds on malignant_markers for a refined classification

# Save the annotated AnnData object with malignant/non-malignant classification
sc.write('./output/cell_annotation/malignant_classification.h5ad', adata)
