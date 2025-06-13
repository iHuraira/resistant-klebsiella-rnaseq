rule bbmap_rRNA:
    input:
        forward_trim_file = "results/trim/{sample}/{sample}_1.fastq",
        reverse_trim_file = "results/trim/{sample}/{sample}_2.fastq"
    conda:
        "../envs/bbmap.yaml"
    output:
        forward_output_file = "results/rRNA/{sample}/{sample}_1.fastq",
        reverse_output_file = "results/rRNA/{sample}/{sample}_2.fastq"
    params:
        k = 31,
        hdist = 1,
        stats = "results/stats/{sample}_bbduk_stats.txt",
        rna_db_path = config["ref_database_rRNA"]
    shell:
        """
        bbduk.sh \
        in1={input.forward_trim_file} \
        in2={input.reverse_trim_file} \
        out1={output.forward_output_file} \
        out2={output.reverse_output_file} \
        ref={params.rna_db_path} \
        k={params.k} \
        hdist={params.hdist} \
        stats={params.stats}
        """
