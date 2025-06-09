rule raw_quality_check:
    input:
        forward_file = "results/fastq/{sample}/{sample}_1.fastq",
        reverse_file = "results/fastq/{sample}/{sample}_2.fastq"
    conda:
        "../envs/quality_control.yaml"
    output:
        forward_html = "results/qc/{sample}/{sample}_1_fastqc.html",
        reverse_html = "results/qc/{sample}/{sample}_2_fastqc.html",
        forward_zip = "results/qc/{sample}/{sample}_1_fastqc.zip",
        reverse_zip = "results/qc/{sample}/{sample}_2_fastqc.zip"
    shell:
        """
        fastqc {input.forward_file} --outdir results/qc/{wildcards.sample}
        fastqc {input.reverse_file} --outdir results/qc/{wildcards.sample}
        """
