#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/001
conda install -y bwa samtools
bwa index ./data/WGS_WES/minigenome.fa
bwa mem -M ./data/WGS_WES/minigenome.fa ./output/001/clean_rnaseq_1.fastq.gz ./output/001/clean_rnaseq_2.fastq.gz > ./output/001/aligned.sam
