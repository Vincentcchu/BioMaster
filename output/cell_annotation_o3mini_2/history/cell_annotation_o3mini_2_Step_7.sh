#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_2
conda install -y -c conda-forge scanpy celltypist pandas
cat << 'EOF' > ./output/cell_annotation_2/annotate_celltypes.py
#!/usr/bin/env python
import scanpy as sc
import pandas as pd
import celltypist
adata = sc.read_h5ad('./output/cell_annotation_2/clustered_adata.h5ad')
markers = pd.read_csv('./output/cell_annotation_2/deg_results.csv', header=0)
celltypist.models.download_models()
predictions = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True)
adata.obs["cell_type"] = predictions.predicted_labels
def classify(cell_type):
    if "Malignant" in cell_type:
        return "malignant"
    else:
        return "non-malignant"
adata.obs["malignancy"] = adata.obs["cell_type"].apply(classify)
adata.write('./output/cell_annotation_2/annotated_adata.h5ad')
EOF
python ./output/cell_annotation_2/annotate_celltypes.py
