rule download_sra:
    output: "results/sra/{sample}/{sample}.sra"
    conda: "../envs/sra_tools.yaml"
    shell:
        """
        prefetch {wildcards.sample} --output-file {output}
        """
