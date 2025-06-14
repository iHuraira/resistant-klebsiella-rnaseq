configfile: "config/config.yaml"

with open(config["sra_path"], "r") as sra_text_file:
    samples = [line.strip() for line in sra_text_file if line.strip()]

include: "rules/sra_tools.smk"
include: "rules/fastq_convert.smk"
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
include: "rules/trim.smk"
include: "rules/post_trim_qc.smk"
include: "rules/post_trim_multiqc.smk"
include: "rules/bbmap.smk"
include: "rules/star_index.smk"
include: "rules/star_align.smk"
include: "rules/RSeQC.smk"
include: "rules/count_reads.smk"
include: "rules/combine_counts.smk"
include: "rules/differential_expression.smk"

rule all:
    input:
        expand("results/sra/{sample}/{sample}.sra", sample=samples),

        expand("results/fastq/{sample}/{sample}_1.fastq", sample=samples),
        expand("results/fastq/{sample}/{sample}_2.fastq", sample=samples),

        expand("results/qc/{sample}/{sample}_1_fastqc.html", sample=samples),
        expand("results/qc/{sample}/{sample}_1_fastqc.zip", sample=samples),
        expand("results/qc/{sample}/{sample}_2_fastqc.html", sample=samples),
        expand("results/qc/{sample}/{sample}_2_fastqc.zip", sample=samples),

        "results/multiqc/multiqc_report.html",

        expand("results/trim/{sample}/{sample}_1.fastq", sample=samples),
        expand("results/trim/{sample}/{sample}_2.fastq", sample=samples),
        expand("results/trim/{sample}/{sample}_1_unpaired.fastq", sample=samples),
        expand("results/trim/{sample}/{sample}_2_unpaired.fastq", sample=samples),

        expand("results/post_trim_qc/{sample}/{sample}_1_fastqc.html", sample=samples),
        expand("results/post_trim_qc/{sample}/{sample}_1_fastqc.zip", sample=samples),
        expand("results/post_trim_qc/{sample}/{sample}_2_fastqc.html", sample=samples),
        expand("results/post_trim_qc/{sample}/{sample}_2_fastqc.zip", sample=samples),

        "results/post_trim_multiqc/multiqc_report.html",

        expand("results/rRNA/{sample}/{sample}_1.fastq", sample=samples),
        expand("results/rRNA/{sample}/{sample}_2.fastq", sample=samples),

        "resources/star_index",

        expand("results/star/{sample}/{sample}.bam", sample=samples),

        expand("results/strandedness/{sample}/{sample}_strandedness.txt", sample=samples),

        expand("results/counts/{sample}/{sample}_counts.txt", sample=samples),

        "results/counts/all_counts.tsv",

        "results/deseq2/deseq2_results.csv",
        "results/deseq2/ma_plot.pdf",
        "results/deseq2/pca_plot.pdf",
        "results/deseq2/sample_heatmap.pdf",
        "results/deseq2/volcano_plot.pdf"

