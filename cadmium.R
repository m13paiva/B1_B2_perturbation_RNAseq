library(DESeq2)
library(tximport)
library(readr)
library(dplyr)
library(ggplot2)

# paths
samples_file <- "/home/m13paiva/Desktop/rna_seq_tese/cadmium_treatment/samples_filtered.tsv"
txi_file <- "/home/m13paiva/Desktop/rna_seq_tese/cadmium_treatment/PRJNA890021/tximport/txi.rds"
targets_file <- "/home/m13paiva/Desktop/rna_seq_tese/cadmium_treatment/targets.tsv"
out_dir <- "/home/m13paiva/Desktop/rna_seq_tese/cadmium_treatment/plots"
out_tsv <- "/home/m13paiva/Desktop/rna_seq_tese/cadmium_treatment/riboflavin_targets_deseq2_results.tsv"

# load samples and txi
samples <- as.data.frame(read_tsv(samples_file, show_col_types = FALSE))
rownames(samples) <- samples$SRR
samples$treatment <- factor(samples$treatment, levels = c("0mg/L", "4mg/L"))
samples$replicate <- factor(samples$replicate)

txi <- readRDS(txi_file)
txi$counts <- txi$counts[, rownames(samples)]
txi$abundance <- txi$abundance[, rownames(samples)]
txi$length <- txi$length[, rownames(samples)]

# deseq2 GLM (~ treatment) with empirical Bayes dispersion shrinkage
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ treatment)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = FALSE)

# 1. exploratory pca plot (top 500 variable genes)
pca_data <- plotPCA(vsd, intgroup = c("treatment", "replicate"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = treatment, fill = treatment, shape = replicate), size = 4.8, stroke = 1.1, alpha = 0.95) +
  scale_color_manual(values = c("0mg/L" = "#2B83BA", "4mg/L" = "#D7191C")) +
  scale_fill_manual(values = c("0mg/L" = "#2B83BA", "4mg/L" = "#D7191C")) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(face = "bold", size = 9.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = sprintf("PC1: %d%% variance", percent_var[1]),
    y = sprintf("PC2: %d%% variance", percent_var[2]),
    color = "Dose (Cd)",
    fill = "Dose (Cd)",
    shape = "Replicate"
  )

ggsave(file.path(out_dir, "cadmium_wt_pca_plot.png"), plot = p_pca, width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "cadmium_wt_pca_plot.pdf"), plot = p_pca, width = 7.5, height = 5.5)

# 2. statistical significance testing: pairwise Wald test (4mg/L vs 0mg/L) + BH FDR adjustment
res <- results(dds, contrast = c("treatment", "4mg/L", "0mg/L"))
targets_all <- read_tsv(targets_file, show_col_types = FALSE)
ribo_targets <- targets_all %>% filter(pathway == "Riboflavin Biosynthesis")

srr_0mg <- samples$SRR[samples$treatment == "0mg/L"]
srr_4mg <- samples$SRR[samples$treatment == "4mg/L"]

group_means <- data.frame(
  locus_id = rownames(norm_counts),
  mean_normalized_0mgL = rowMeans(norm_counts[, srr_0mg, drop = FALSE]),
  mean_normalized_4mgL = rowMeans(norm_counts[, srr_4mg, drop = FALSE]),
  stringsAsFactors = FALSE
)

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column(var = "locus_id") %>%
  rename(
    log2FC_4vs0 = log2FoldChange,
    lfcSE_4vs0 = lfcSE,
    stat_4vs0 = stat,
    pvalue_4vs0 = pvalue,
    padj_4vs0 = padj
  )

ribo_res <- ribo_targets %>%
  left_join(res_df, by = "locus_id") %>%
  left_join(group_means, by = "locus_id") %>%
  mutate(
    log2FoldChange = log2FC_4vs0,
    lfcSE = lfcSE_4vs0,
    pvalue = pvalue_4vs0,
    padj = padj_4vs0,
    passed_filter = !is.na(baseMean) & baseMean >= 10,
    filter_status = ifelse(passed_filter, "PASS", "FILTERED_LOW_EXPRESSION"),
    expression_status = ifelse(passed_filter, "Adequate Expression (baseMean >= 10)", "Low/Unexpressed (baseMean < 10)"),
    significance_flag = case_when(
      !passed_filter ~ "Low Expression",
      is.na(padj) ~ "ns",
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    significance_flag_4vs0 = significance_flag,
    regulation_category = factor(
      case_when(
        !passed_filter ~ "Low Expression (baseMean < 10)",
        is.na(padj) | padj >= 0.05 ~ "Not Significant (padj >= 0.05)",
        log2FoldChange > 0 ~ "Upregulated (padj < 0.05)",
        log2FoldChange < 0 ~ "Downregulated (padj < 0.05)"
      ),
      levels = c(
        "Upregulated (padj < 0.05)",
        "Downregulated (padj < 0.05)",
        "Not Significant (padj >= 0.05)",
        "Low Expression (baseMean < 10)"
      )
    ),
    plot_label = identifier,
    plot_label = factor(plot_label, levels = plot_label[order(log2FoldChange, na.last = FALSE)]),
    star_pos = ifelse(
      log2FoldChange >= 0,
      log2FoldChange + ifelse(!is.na(lfcSE), lfcSE, 0) + 0.20,
      log2FoldChange - ifelse(!is.na(lfcSE), lfcSE, 0) - 0.20
    )
  )

write_tsv(ribo_res, out_tsv)

# 3. riboflavin log2fc summary barplot
regulation_colors <- c(
  "Upregulated (padj < 0.05)" = "#D73027",
  "Downregulated (padj < 0.05)" = "#4575B4",
  "Not Significant (padj >= 0.05)" = "#999999",
  "Low Expression (baseMean < 10)" = "#E0E0E0"
)

p_bar <- ggplot(ribo_res, aes(x = plot_label, y = log2FoldChange, fill = regulation_category)) +
  geom_col(color = "black", linewidth = 0.35, width = 0.75) +
  geom_errorbar(
    data = ribo_res %>% filter(passed_filter == TRUE),
    aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
    width = 0.35, linewidth = 0.45, color = "gray30"
  ) +
  geom_text(
    data = ribo_res %>% filter(passed_filter == TRUE & significance_flag != "ns"),
    aes(y = star_pos, label = significance_flag),
    size = 4.2, fontface = "bold", color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = regulation_colors, drop = FALSE) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "bold", size = 8.5),
    axis.text.x = element_text(face = "bold", size = 8.5),
    axis.title = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.5),
    panel.grid.major.y = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  ) +
  labs(
    x = "Target Gene",
    y = "Log2 Fold Change (4mg/L Cd / 0mg/L Cd)",
    fill = "Differential Expression Status"
  )

ggsave(file.path(out_dir, "riboflavin_log2fc_summary_barplot.png"), plot = p_bar, width = 5.5, height = 9.5, dpi = 300)
ggsave(file.path(out_dir, "riboflavin_log2fc_summary_barplot.pdf"), plot = p_bar, width = 5.5, height = 9.5)
