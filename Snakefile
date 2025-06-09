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
