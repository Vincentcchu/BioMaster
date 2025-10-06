#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation
conda install -y scanpy celltypist pandas
python ./output/cell_annotation/malignant_classification.py --dataset ./output/cell_annotation/normalized_dataset.h5ad --clusters ./output/cell_annotation/clusters_annotated.h5ad --markers ./output/cell_annotation/cluster_markers.csv --output ./output/cell_annotation/malignant_classification.h5ad
