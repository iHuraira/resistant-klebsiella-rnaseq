rule star_index:
    input:
        genome_file = config["genome_path"],
        annotations_file = config["annotations_path"]
    conda:
        "../envs/star.yaml"
    output:
        output_dir = directory("resources/star_index")
    params:
        threads = 6
    shell:
        """
        STAR \
        --runThreadN {params.threads} \
        --runMode genomeGenerate \
        --genomeDir {output.output_dir} \
        --genomeFastaFiles {input.genome_file} \
        --sjdbGTFfile {input.annotations_file} \
        --sjdbOverhang 99
        """
