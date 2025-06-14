rule combine_counts:
    input:
        expand("results/counts/{sample}/{sample}_counts.txt", sample=samples)
    output:
        "results/counts/all_counts.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_counts.py"
