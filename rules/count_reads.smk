rule count_reads:
    input:
        bam = "results/star/{sample}/{sample}.bam",
        gtf = "resources/Genome/annotations.gtf"
    output:
        counts = "results/counts/{sample}/{sample}_counts.txt"
    conda:
        "../envs/subread.yaml"
    params:
        threads = 4,
        stranded = 0  # 0 = unstranded
    shell:
        """
        featureCounts \
            -T {params.threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -s {params.stranded} \
            -p \
            {input.bam}
        """
