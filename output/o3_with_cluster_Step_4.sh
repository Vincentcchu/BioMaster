#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_with_cluster
conda install -y scanpy anndata pandas scikit-learn
mkdir -p ./output/o3_with_cluster
grep -v "^$" ./output/o3_with_cluster/marker_genes.txt > ./output/o3_with_cluster/marker_genes_clean.txt
if [ ! -s ./output/o3_with_cluster/marker_genes_clean.txt ]; then echo Marker file is empty; exit 1; fi
python ./output/o3_with_cluster/classify_cells.py --input ./output/o3_with_cluster/clustered_dataset.h5ad --markers ./output/o3_with_cluster/marker_genes_clean.txt --output ./output/o3_with_cluster/classified_cells.h5ad
mv ./output/o3_with_cluster/classified_cells.h5ad ./output/o3_with_cluster/cell_classification_results.h5ad
python ./output/o3_with_cluster/export_classification.py --input ./output/o3_with_cluster/cell_classification_results.h5ad --output ./output/o3_with_cluster/classified_cells.txt
