rule trimming:
    input:
        forward_file = "results/fastq/{sample}/{sample}_1.fastq",
        reverse_file = "results/fastq/{sample}/{sample}_2.fastq"
    conda:
        "../envs/trim.yaml"
    params:
        adapters = config["adapters_path"]
    output:
        forward_trim = "results/trim/{sample}/{sample}_1.fastq",
        reverse_trim = "results/trim/{sample}/{sample}_2.fastq",
        forward_unpaired = "results/trim/{sample}/{sample}_1_unpaired.fastq",
        reverse_unpaired = "results/trim/{sample}/{sample}_2_unpaired.fastq"
    threads: 6
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 \
        {input.forward_file} {input.reverse_file} \
        {output.forward_trim} {output.forward_unpaired} \
        {output.reverse_trim} {output.reverse_unpaired} \
        ILLUMINACLIP:{params.adapters}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """