#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/001
conda install -y -c bioconda gatk4 samtools
samtools faidx ./data/WGS_WES/minigenome.fa
gatk CreateSequenceDictionary -R ./data/WGS_WES/minigenome.fa -O ./data/WGS_WES/minigenome.dict
gatk BaseRecalibrator -R ./data/WGS_WES/minigenome.fa -I ./output/001/aligned.dedup.bam --known-sites ./data/WGS_WES/known_sites.vcf -O ./output/001/recal_data.table
gatk ApplyBQSR -R ./data/WGS_WES/minigenome.fa -I ./output/001/aligned.dedup.bam --bqsr-recal-file ./output/001/recal_data.table -O ./output/001/aligned.dedup.recal.bam
samtools index ./output/001/aligned.dedup.recal.bam
