#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/001
conda install -y gatk4 samtools openjdk
gatk MarkDuplicates \  -I ./output/001/aligned.sorted.bam \  -O ./output/001/aligned.dedup.bam \  -M ./output/001/marked_dup_metrics.txt \  --REMOVE_DUPLICATES false \  --VALIDATION_STRINGENCY SILENT
samtools index ./output/001/aligned.dedup.bam
