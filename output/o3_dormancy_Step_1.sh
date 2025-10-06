#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_dormancy
conda install -y scanpy
mkdir -p ./output/o3_dormancy
cat << 'EOF' > ./output/o3_dormancy/subset_malignant.py
import scanpy as sc

def main():
    # Load the malignant single-cell RNA-seq dataset
    input_file = '/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/biomaster/output/cell_annotation_o3mini/malignant_classification.h5ad'
    adata = sc.read_h5ad(input_file)
    
    # Save the AnnData object to the designated output file
    output_file = './output/o3_dormancy/initialized_malignant_cells.h5ad'
    adata.write(output_file)

if __name__ == '__main__':
    main()
EOF
python ./output/o3_dormancy/subset_malignant.py
