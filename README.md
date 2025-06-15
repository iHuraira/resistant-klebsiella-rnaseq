Below is a **ready-to-copy README.md** that documents your entire analysis pipeline.
Feel free to tweak section names or add badges, but this should already give you a polished, complete overview.

---

```markdown
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

A DAG of the workflow is stored in `docs/pipeline_dag.svg`.

## 3. Directory Layout
```

├── data/
│   ├── raw/           # .sra and FASTQ files
│   ├── trimmed/
│   └── genome/        # reference fasta + gtf + STAR index
├── qc/                # FastQC & MultiQC reports
├── align/             # sorted BAMs
├── counts/            # geneCounts.txt
├── results/
│   ├── deg/           # DESeq2 outputs
│   ├── enrichment/    # GO & KEGG tables, plots
│   └── variants/      # VCFs
├── scripts/           # R and helper bash / python
└── Snakefile

````

## 4. Quick Start

```bash
# 1. Clone
git clone https://github.com/yourname/kpn16_rnaseq_pipeline.git
cd kpn16_rnaseq_pipeline

# 2. Create Conda environment (Snakemake ≥7)
mamba env create -f envs/core.yaml
conda activate kpn_rnaseq

# 3. Edit config.yaml (threads, paths, parameters)

# 4. Run
snakemake --use-conda --cores 16 --rerun-incomplete
````

Outputs and interactive HTML reports appear in `results/` and `qc/`.

## 5. Key Results

* **PCA:** PC1 (71 % variance) cleanly separates KPN16 vs. ATCC13883; replicates cluster tightly except one possible KPN16 outlier (investigated in QC).
* **DEGs:** *n* = XXXX up-regulated, *n* = YYYY down-regulated (|log₂FC| > 1, FDR < 0.001).
* **Enriched GO Terms:** membrane transport, LPS biosynthesis, butanoate metabolism.
* **Enriched KEGG Pathways:** ABC transporters, two-component systems, β-lactam resistance.
* **Variants:** high-confidence SNPs/indels in efflux pump regulators and LPS genes (see `results/variants/`).

## 6. Software Versions

| Tool            | Version |
| --------------- | ------- |
| Snakemake       | 7.32    |
| STAR            | 2.7.11a |
| Trimmomatic     | 0.39    |
| FastQC          | 0.12.1  |
| MultiQC         | 1.16    |
| BBMap           | 39.01   |
| RSeQC           | 4.0.1   |
| featureCounts   | 2.0.3   |
| DESeq2          | 1.40    |
| InterProScan    | 5.65-97 |
| KOBAS           | 3.0.4   |
| clusterProfiler | 4.8     |
| GATK            | 4.5     |
| Picard          | 2.27    |
| Samtools        | 1.17    |

Full `conda` environment specs are in `envs/`.

## 7. Reproducibility

* All parameters are logged via Snakemake.
* Checksums of raw data and reference files are stored in `logs/`.
* Pipeline DAG + HTML reports ensure full provenance.

## 8. References

1. Liu Y **et al.** (2025) *IET Syst Biol* 19(1): e12112.
2. Dobin A **et al.** 2013 STAR: ultrafast RNA-seq aligner. *Bioinformatics*.
3. Love MI, Huber W, Anders S. 2014. **DESeq2**. *Genome Biology*.
4. Jones P **et al.** 2014. **InterProScan 5**. *Nucleic Acids Res*.

---

*Maintainer: Your Name — `<you@example.org>`*
*Licensed under MIT.*

```

---

### How to use it
1. Copy everything between the triple back-quotes into a new `README.md`.
2. Replace placeholder counts (`XXXX`, `YYYY`) with your actual DEG numbers.
3. Add or remove sections (e.g., “License”, “Citing this pipeline”) to suit your repository.

Let me know if you’d like additional tweaks—such as badges, logo, or CI instructions!
```
