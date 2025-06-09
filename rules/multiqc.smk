rule multi_quality_check:
    input:
        expand("results/qc/{sample}/{sample}_1_fastqc.zip", sample=samples),
        expand("results/qc/{sample}/{sample}_2_fastqc.zip", sample=samples)
    conda:
        "../envs/quality_control.yaml"
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        """
        multiqc results/qc/ -o ./results/multiqc
        """