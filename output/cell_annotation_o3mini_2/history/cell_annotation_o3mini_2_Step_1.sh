#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_2
conda install -y scanpy
python -c "import scanpy as sc; adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad'); adata.write('./output/cell_annotation_2/initialized_dataset.h5ad')"
