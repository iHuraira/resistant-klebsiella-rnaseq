# 🧬 Snakemake Pipeline for RNA-Seq Analysis of *Klebsiella pneumoniae* Antibiotic Resistance

## Introduction

This repository contains a Snakemake pipeline designed for RNA-Seq analysis of *Klebsiella pneumoniae* isolates to study mechanisms of antibiotic resistance. The analysis is based on the publicly available dataset [**GSE229867**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229867), which profiles gene expression in two bacterial strains:

* **ATCC13883** – reference strain
* **KPN16** – multidrug-resistant clinical isolate

The original study, published in *IET Systems Biology* (Liu et al., 2025), revealed that resistance in KPN16 may be associated with increased butanoate metabolism and lipopolysaccharide biosynthesis, alongside reduced transmembrane transport activity.

The pipeline enables automated, reproducible processing of RNA-Seq data including quality control, alignment, quantification, and differential expression analysis.

Great! Here's a well-formatted **Pipeline Overview** section for your README that clearly describes each step while maintaining clarity and conciseness:

---

## 🧪 Pipeline Overview

This Snakemake pipeline automates the analysis of RNA-Seq data from *Klebsiella pneumoniae*, covering raw data retrieval to differential expression analysis. The steps are modular and reproducible, suitable for scaling and customization.

### Workflow Steps:

1. **prefetch\_sra**
   Downloads raw sequencing data from the NCBI Sequence Read Archive (SRA) using `prefetch`.

2. **convert\_to\_fastq**
   Converts `.sra` files into `.fastq` format using `fasterq-dump`.

3. **trim\_reads**
   Trims adapter sequences and low-quality bases using **Trimmomatic**.

4. **qc\_prepost**
   Performs quality control on raw and trimmed reads using **FastQC**, then compiles reports with **MultiQC**.

5. **check\_rRNA\_content**
   Estimates ribosomal RNA contamination levels using **BBMap** (`bbduk.sh`).

6. **align\_reads**
   Aligns reads to the reference genome using **STAR** aligner in 2-pass mode for improved accuracy.

7. **check\_strandedness**
   Assesses library strandedness using **RSeQC** tools (e.g., `infer_experiment.py`).

8. **count\_reads**
   Counts aligned reads at the gene level using **featureCounts**.

9. **differential\_expression**
   Identifies differentially expressed genes between resistant and non-resistant strains using **DESeq2** or **edgeR** in R.

10. **Gene Set Enrichment Analysis (GSEA)** *(performed independently)*
    Functional enrichment analysis was conducted using **InterProScan**, although it is not part of the automated pipeline.

---
