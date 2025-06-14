library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(stringr)

# === Load Data ===
counts <- read_tsv(snakemake@input[["counts"]])
metadata <- read.csv(snakemake@input[["metadata"]], row.names = 1)

# Sanitize metadata and counts
metadata$condition <- factor(metadata$condition)
colnames(counts)[-1] <- str_trim(colnames(counts)[-1])
rownames(metadata) <- str_trim(rownames(metadata))

# Prepare count matrix
count_matrix <- as.matrix(counts[,-1])
rownames(count_matrix) <- counts[[1]]

print("Count matrix column names:")
print(colnames(count_matrix))

print("Metadata rownames:")
print(rownames(metadata))

print("Metadata sample column (if available):")
print(metadata$sample)

# Ensure sample names match
stopifnot(all(colnames(count_matrix) %in% rownames(metadata)))

# Match metadata order
metadata <- metadata[colnames(count_matrix), , drop = FALSE]

# === DESeq2 Analysis ===
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res) %>%
  mutate(gene = rownames(.)) %>%
  arrange(padj)

# Save results
write.csv(res_df, snakemake@output[["results"]], row.names = FALSE)

# === MA Plot ===
pdf(snakemake@output[["ma"]])
plotMA(res, ylim = c(-5, 5))
dev.off()

# === PCA Plot ===
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_minimal()
ggsave(snakemake@output[["pca"]], p, width = 6, height = 5)

# === Heatmap of Top Genes ===
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

top_n <- min(30, nrow(top_genes))
top_genes <- top_genes %>%
  slice_head(n = top_n) %>%
  pull(gene)

mat <- assay(vsd)[top_genes, , drop = FALSE]
mat <- mat - rowMeans(mat)
pdf(snakemake@output[["heatmap"]])
pheatmap(mat,
         annotation_col = metadata,
         show_rownames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))
dev.off()

# === Volcano Plot ===
res_df <- res_df %>% filter(!is.na(padj))
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "yes", "no")
res_df$log10padj <- -log10(res_df$padj + 1e-10)
v <- ggplot(res_df, aes(log2FoldChange, log10padj)) +
  geom_point(aes(color = significant), alpha = 0.5) +
  scale_color_manual(values = c("no" = "gray", "yes" = "red")) +
  theme_minimal() +
  ggtitle("Volcano Plot") +
  xlab("log2 Fold Change") + ylab("-log10(padj)")
ggsave(snakemake@output[["volcano"]], v, width = 6, height = 5)
