#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_sonnet3.7t
conda install -c conda-forge -c bioconda -y scanpy matplotlib
mkdir -p ./output/cell_annotation_sonnet3.7t
python -c "import scanpy as sc; import os; os.makedirs('./output/cell_annotation_sonnet3.7t', exist_ok=True); adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad'); adata.var['mt'] = adata.var_names.str.startswith('mt-'); sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True); adata = adata[adata.obs.pct_counts_mt <= 20, :]; adata = adata[adata.obs.total_counts >= 500, :]; sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata); sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5); adata.write('./output/cell_annotation_sonnet3.7t/preprocessed_data.h5ad'); print('Preprocessing complete.')"
