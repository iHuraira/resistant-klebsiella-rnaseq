# Klebsiella pneumoniae RNA-seq Pipeline  
*Antibiotic-resistance transcriptomics (GSE229867)*

## 1. Project Overview
The increasing antibiotic resistance of **_Klebsiella pneumoniae_** poses a serious threat to global public health.  
To investigate resistance mechanisms, we re-analysed the RNA-seq data from GEO Series **GSE229867**, comparing the clinical mutant strain **KPN16** with the reference strain **ATCC13883** (three biological replicates each).

> **Original study**  
> Liu Y *et al.* 2025. *Transcriptomic analysis reveals pathways underlying the multi-antibiotic resistance of Klebsiella pneumoniae.* **IET Syst Biol** 19(1): e12112.  

### Experimental design (from GEO)
- **Platform:** Illumina NovaSeq 6000 (paired-end 150 bp)  
- **BioProject:** PRJNA956495  
- **Samples:**  
  | Condition | GEO IDs | Replicates |
  |-----------|---------|------------|
  | ATCC13883 | GSM7179288-90 | 3 |
  | KPN16     | GSM7179288-90 | 3 |

Each isolate was cultured to OD<sub>600</sub>=0.05, total RNA extracted, rRNA removed, libraries built (370–420 bp inserts) and sequenced on NovaSeq.

## 2. Pipeline Summary
Implemented in **Snakemake** (`Snakefile`) with modular rules:

| Order | Rule / Tool | Purpose |
|-------|-------------|---------|
| 1 | `prefetch_sra` (– NCBI SRA) | Fetch .sra files |
| 2 | `convert_to_fastq` (`fasterq-dump`) | Convert to paired-end FASTQ |
| 3 | `trim_reads` (**Trimmomatic**) | Adapter & low-quality trimming |
| 4 | `qc_prepost` (**FastQC**, **MultiQC**) | Pre- & post-trim quality reports |
| 5 | `check_rRNA_content` (**BBMap**) | rRNA contamination screening |
| 6 | `align_reads` (**STAR**) | Genome alignment to GCA_000240185.2 |
| 7 | `check_strandedness` (**RSeQC**) | Library strand orientation |
| 8 | `count_reads` (**featureCounts**) | Gene-level read counts |
| 9 | `differential_expression` (**DESeq2**) | Identify DEGs (|log₂FC| ≥ 1, FDR ≤ 0.001) |
| 10 | `interproscan_annotations` (**InterProScan 5**) | Assign GO terms |
| 11 | `kegganalysis` (**KOBAS 3**) | KEGG pathway mapping & enrichment |
| 12 | `go_kegg_enrichment` (**clusterProfiler**) | GO/KEGG over-representation tests |
| 13 | `variant_calling` (**GATK HaplotypeCaller → CombineGVCFs → GenotypeGVCFs**, filtered with MQ < 40, QD < 2, FS > 30, DP < 10, QUAL < 20) |
| 14 | `visualisation` (R ggplot2, pheatmap) | PCA, volcano plots, heatmaps |


## 8. References

1. Liu Y **et al.** (2025) *IET Syst Biol* 19(1): e12112.
2. Dobin A **et al.** 2013 STAR: ultrafast RNA-seq aligner. *Bioinformatics*.
3. Love MI, Huber W, Anders S. 2014. **DESeq2**. *Genome Biology*.
4. Jones P **et al.** 2014. **InterProScan 5**. *Nucleic Acids Res*.
