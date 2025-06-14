rule check_strandedness:
    input:
        bam = "results/star/{sample}/{sample}.bam",
        bed = "resources/annotations.bed"
    output:
        txt = "results/strandedness/{sample}/{sample}_strandedness.txt"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        infer_experiment.py \
            -r {input.bed} \
            -i {input.bam} \
            > {output.txt}
        """
