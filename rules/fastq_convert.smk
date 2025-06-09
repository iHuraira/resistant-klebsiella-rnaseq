rule fastq_files_converter:
    input:
        "results/sra/{sample}/{sample}.sra"
    output:
        forward_file = "results/fastq/{sample}/{sample}_1.fastq",
        reverse_file = "results/fastq/{sample}/{sample}_2.fastq"
    conda:
        "../envs/sra_tools.yaml"
    threads: 6
    shell:
        """
        fasterq-dump {wildcards.sample} -e {threads} -O results/fastq/{wildcards.sample}/
        """
