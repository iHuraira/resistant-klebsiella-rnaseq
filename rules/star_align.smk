rule star_align:
    input:
        forward_file = "results/rRNA/{sample}/{sample}_1.fastq",
        reverse_file = "results/rRNA/{sample}/{sample}_2.fastq",
        index_dir = "resources/star_index"
    output:
        bam = "results/star/{sample}/{sample}.bam",
        log = "results/star/{sample}/{sample}_Log.final.out",
        counts = "results/star/{sample}/{sample}_ReadsPerGene.out.tab"
    params:
        threads = 3,
        out_prefix = "results/star/{sample}/{sample}_"
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR \
            --runThreadN {params.threads} \
            --genomeDir {input.index_dir} \
            --readFilesIn {input.forward_file} {input.reverse_file} \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM \
            --outSAMunmapped Within \
            --twopassMode Basic \
            --quantMode GeneCounts

        # Rename only the BAM file
        mv {params.out_prefix}Aligned.sortedByCoord.out.bam {output.bam}

        # Clean up intermediate/unnecessary files
        rm -f {params.out_prefix}Log.out \
              {params.out_prefix}Log.progress.out \
              {params.out_prefix}SJ.out.tab

        rm -rf {params.out_prefix}__STARpass1 \
               {params.out_prefix}__STARgenome
        """
