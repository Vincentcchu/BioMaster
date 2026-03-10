# Workflow for Single-Cell RNA-sequencing (scRNA-seq) Experimental and Computational Pipeline

The following step-by-step workflow summarizes the typical procedures described in Haque et al. (2017) for designing, executing, and analyzing an scRNA-seq study. Each step focuses on one tool or procedure and clearly states its purpose, inputs, and expected outputs.

---

## Step 1: Single-Cell Isolation  
**Description:** Isolate viable individual cells from the tissue sample while minimizing stress‐induced transcriptomic changes.  
**Input Required:** Fresh tissue sample or dissociated cells; cell suspension.  
**Expected Output:** A suspension of single cells ready for further processing.  
**Tool Used:** Fluorescence Activated Cell Sorting (FACS)  
*(Note: FACS is an example method; alternative micro-dissection or microfluidic platforms may be used.)*

---

## Step 2: Cell Viability Assessment  
**Description:** Assess the quality and viability of isolated cells to ensure that only live, intact cells proceed to downstream processing.  
**Input Required:** Isolated single cell suspension; viability dyes (e.g., propidium iodide, Calcein AM).  
**Expected Output:** A record/report of viable versus nonviable cells, with a purified pool of viable single cells.  
**Tool Used:** Flow Cytometer (integrated with FACS or as a stand-alone viability assay)

---

## Step 3: Cell Lysis  
**Description:** Lyse each single cell under conditions that preserve mRNA integrity for subsequent capture.  
**Input Required:** Collection of isolated, viable single cells; appropriate cell lysis buffer/reagents.  
**Expected Output:** Lysed cells with released, intact mRNA molecules.  
**Tool Used:** Cell Lysis Kit (optimized for single-cell applications)

---

## Step 4: mRNA Capture  
**Description:** Capture polyadenylated mRNA molecules from the lysate using oligo[dT] primers to selectively target mRNA over ribosomal RNA.  
**Input Required:** Lysed cell solution containing mRNA; magnetic beads or plate-bound oligo[dT] primers.  
**Expected Output:** mRNA molecules bound to capture reagents, ready for reverse transcription.  
**Tool Used:** Poly(T) Bead-Based Capture System

---

## Step 5: Reverse Transcription  
**Description:** Convert the captured mRNA into complementary DNA (cDNA) using reverse transcriptase, with the addition of adaptor sequences and optionally Unique Molecular Identifiers (UMIs) to tag each transcript.  
**Input Required:** Captured mRNA; reverse transcription reagents (enzymes, primers, dNTPs).  
**Expected Output:** First-strand cDNA with incorporated adaptor and UMI sequences.  
**Tool Used:** SMARTer Reverse Transcriptase (or a similar reverse transcription kit optimized for scRNA-seq)

---

## Step 6: cDNA Amplification  
**Description:** Amplify the minute quantities of cDNA obtained from individual cells using PCR to generate sufficient material for library preparation.  
**Input Required:** First-strand cDNA; PCR reagents (polymerase, primers specific for adaptor sequences).  
**Expected Output:** An amplified cDNA pool enriched in transcripts from each cell.  
**Tool Used:** Thermal Cycler (PCR Machine)

---

## Step 7: Library Preparation (Barcode-tagging and Adapter Ligation)  
**Description:** Prepare sequencing libraries by ligating sequencing adapters and incorporating cellular barcodes that allow sample multiplexing and assignment of reads to individual cells.  
**Input Required:** Amplified cDNA; library preparation reagents and barcoding/adaptor kits.  
**Expected Output:** Barcoded cDNA libraries ready for next-generation sequencing (NGS).  
**Tool Used:** Commercial Library Preparation Kit (e.g., Illumina Nextera XT or a scRNA-seq–specific barcoding kit)

---

## Step 8: Pooling and Next-Generation Sequencing  
**Description:** Combine barcoded libraries from multiple cells and sequence the pooled library on an NGS platform.  
**Input Required:** Individual barcoded cDNA libraries; pooled library preparation.  
**Expected Output:** Raw sequencing reads (FASTQ files) containing cell-specific barcode and UMI information.  
**Tool Used:** Illumina Sequencing Platform (or equivalent NGS instrument)

---

## Step 9: Quality Control of Raw Sequencing Data  
**Description:** Evaluate the quality of raw sequencing data to detect issues such as low quality scores, adapter contamination, or low read complexity.  
**Input Required:** Raw FASTQ sequencing files.  
**Expected Output:** Quality control reports (e.g., per-base quality metrics, GC content, duplication levels).  
**Tool Used:** FastQC

---

## Step 10: Read Alignment to Reference Genome  
**Description:** Map sequencing reads to a reference genome or transcriptome to determine the genomic origin of each read, while preserving cell barcode and UMI information.  
**Input Required:** Raw sequencing reads (FASTQ files); reference genome/transcriptome.  
**Expected Output:** Aligned reads in BAM/SAM file format with cell barcode annotations.  
**Tool Used:** STAR Aligner (or another RNA-seq aligner appropriate for single-cell data)

---

## Step 11: Transcript Quantification and UMI Deduplication  
**Description:** Quantify transcript abundances for individual cells by counting reads/UMIs, reducing amplification bias by collapsing duplicates.  
**Input Required:** Aligned BAM/SAM files; gene annotation files.  
**Expected Output:** Count matrices that indicate the number of transcripts/UMIs per gene for each cell.  
**Tool Used:** FeatureCounts (or UMI-aware quantification tool such as UMI-tools)

---

## Step 12: Data Filtering and Normalization  
**Description:** Filter out poor-quality cells (e.g., with low library size, high mitochondrial gene content) and normalize transcript counts to account for library depth and technical variability.  
**Input Required:** Raw gene count matrix; quality metrics per cell.  
**Expected Output:** A cleaned and normalized gene expression matrix suitable for downstream analyses.  
**Tool Used:** Seurat (or an R-based pipeline such as Scater for pre-processing)

---

## Step 13: Dimensionality Reduction and Clustering  
**Description:** Reduce the high-dimensional gene expression data and perform clustering to identify cell subpopulations, using methods such as PCA and t-SNE.  
**Input Required:** Normalized gene expression matrix.  
**Expected Output:** Low-dimensional representations (e.g., PCA or t-SNE plots) and cluster assignments for each single cell.  
**Tool Used:** Seurat (with built-in PCA/t-SNE modules)

---

## Step 14: Differential Expression and Biological Interpretation  
**Description:** Identify genes that are differentially expressed between clusters or conditions and interpret their biological significance to relate transcriptomic variations to cellular phenotypes.  
**Input Required:** Cluster assignments; normalized gene expression matrix.  
**Expected Output:** Lists of differentially expressed genes, statistical significance values, and enrichment analyses linking gene expression changes to underlying biological processes.  
**Tool Used:** DESeq2 (or Seurat’s differential expression testing functions)

---

This workflow captures the comprehensive steps—from cell isolation to computational downstream analysis—needed for a successful scRNA-seq study as described in Haque et al. (2017). Each individual step is designed to ensure high data quality and reliable biological interpretation, according to the principles outlined in the paper.