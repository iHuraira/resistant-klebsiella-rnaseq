rule post_trim_multi_quality_check:
    input:
        expand("results/post_trim_qc/{sample}/{sample}_1_fastqc.zip", sample=samples),
        expand("results/post_trim_qc/{sample}/{sample}_2_fastqc.zip", sample=samples)
    conda:
        "../envs/quality_control.yaml"
    output:
        "results/post_trim_multiqc/multiqc_report.html"
    shell:
        """
        multiqc results/post_trim_qc/ -o ./results/post_trim_multiqc
        """