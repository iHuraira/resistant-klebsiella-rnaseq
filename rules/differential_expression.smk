rule differential_expression:
    input:
        counts = "results/counts/all_counts.tsv",
        metadata = "resources/samplesheet.csv"
    output:
        results = "results/deseq2/deseq2_results.csv",
        ma = "results/deseq2/ma_plot.pdf",
        pca = "results/deseq2/pca_plot.pdf",
        heatmap = "results/deseq2/sample_heatmap.pdf",
        volcano = "results/deseq2/volcano_plot.pdf"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_analysis.R"
