library(DESeq2)
library(tximport)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(cowplot)

# paths
samples_file <- "/home/m13paiva/Desktop/rna_seq_tese/circadian_rythm/samples_filtered.tsv"
txi_file <- "/home/m13paiva/Desktop/rna_seq_tese/circadian_rythm/PRJNA601442/tximport/txi.rds"
targets_file <- "/home/m13paiva/Desktop/rna_seq_tese/circadian_rythm/targets.tsv"
out_dir <- "/home/m13paiva/Desktop/rna_seq_tese/circadian_rythm/plots"
pub_dir <- "/home/m13paiva/Desktop/rna_seq_tese/circadian_rythm/plots/publication_figures"
thesis_fig_dir <- "/home/m13paiva/Desktop/tese_template/5-Figures/Results"
out_tsv <- "/home/m13paiva/Desktop/rna_seq_tese/circadian_rythm/thiamine_targets_deseq2_results.tsv"

# load samples and txi
samples <- as.data.frame(read_tsv(samples_file, show_col_types = FALSE))
rownames(samples) <- samples$SRR
time_levels <- c("0h", "4h", "8h", "12h", "16h", "20h")
samples$time <- factor(samples$time, levels = time_levels)
samples$time_num <- as.numeric(gsub("h", "", as.character(samples$time)))

txi <- readRDS(txi_file)
txi$counts <- txi$counts[, rownames(samples)]
txi$abundance <- txi$abundance[, rownames(samples)]
txi$length <- txi$length[, rownames(samples)]

# deseq2 GLM (~ time) with empirical Bayes dispersion shrinkage
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ time)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

# 1. exploratory pca plot (top 500 variable genes)
pca_data <- plotPCA(vsd, intgroup = "time", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

circ_colors <- c(
  "0h" = "#2B83BA",
  "4h" = "#74ADD1",
  "8h" = "#FEE090",
  "12h" = "#FDAE61",
  "16h" = "#F46D43",
  "20h" = "#313695"
)

centroids <- pca_data %>%
  group_by(time) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop") %>%
  arrange(match(time, time_levels))
centroids_loop <- bind_rows(centroids, centroids[1, ])

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = time, fill = time)) +
  geom_path(data = centroids_loop, aes(x = PC1, y = PC2), color = "gray55",
            linetype = "dashed", linewidth = 0.85, inherit.aes = FALSE) +
  geom_point(size = 5.2, shape = 21, stroke = 1.15, alpha = 0.95) +
  scale_color_manual(values = circ_colors) +
  scale_fill_manual(values = circ_colors) +
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
    color = "Timepoint",
    fill = "Timepoint"
  )

ggsave(file.path(out_dir, "circadian_timecourse_pca_plot.png"), plot = p_pca, width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "circadian_timecourse_pca_plot.pdf"), plot = p_pca, width = 7.5, height = 5.5)

# 2. statistical significance testing:
# A) Omnibus Likelihood Ratio Test (LRT: ~ time vs ~ 1, df = 5)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)
res_lrt_df <- as.data.frame(res_lrt) %>%
  rownames_to_column(var = "locus_id") %>%
  rename(LRT_stat = stat, pvalue_LRT = pvalue, padj_LRT = padj)

# B) Pairwise Wald tests vs 0h baseline for each timepoint
timepoints_vs_0 <- c("4h", "8h", "12h", "16h", "20h")
res_pairwise <- list()

for (tp in timepoints_vs_0) {
  res_tp <- results(dds, contrast = c("time", tp, "0h"))
  df_tp <- as.data.frame(res_tp) %>%
    rownames_to_column(var = "locus_id") %>%
    select(locus_id, log2FoldChange, lfcSE, pvalue, padj)
  colnames(df_tp)[2:5] <- paste0(c("log2FC_", "lfcSE_", "pvalue_", "padj_"), tp, "_vs_0h")
  res_pairwise[[tp]] <- df_tp
}

all_pairwise_df <- res_pairwise[[1]]
for (k in 2:length(res_pairwise)) {
  all_pairwise_df <- all_pairwise_df %>% left_join(res_pairwise[[k]], by = "locus_id")
}

# C) Mean normalized counts per timepoint and circadian amplitude
targets <- read_tsv(targets_file, show_col_types = FALSE)

time_means_list <- list()
for (tp in time_levels) {
  srr_tp <- samples$SRR[samples$time == tp]
  time_means_list[[paste0("mean_norm_", tp)]] <- rowMeans(norm_counts[, srr_tp, drop = FALSE])
}
time_means_df <- as.data.frame(time_means_list) %>%
  rownames_to_column(var = "locus_id")

mean_mat <- as.matrix(time_means_df[, paste0("mean_norm_", time_levels)])
colnames(mean_mat) <- time_levels
peak_indices <- apply(mean_mat, 1, which.max)
trough_indices <- apply(mean_mat, 1, which.min)

time_summary_df <- data.frame(
  locus_id = time_means_df$locus_id,
  peak_time = time_levels[peak_indices],
  trough_time = time_levels[trough_indices],
  circadian_amplitude_fold = round(ifelse(apply(mean_mat, 1, min) > 0, apply(mean_mat, 1, max) / apply(mean_mat, 1, min), NA), 2),
  stringsAsFactors = FALSE
)

# D) Combine target definitions, LRT statistics, pairwise Wald tests, and baseline filter (baseMean >= 10)
thiamine_res <- targets %>%
  left_join(res_lrt_df %>% select(locus_id, baseMean, LRT_stat, pvalue_LRT, padj_LRT), by = "locus_id") %>%
  left_join(all_pairwise_df, by = "locus_id") %>%
  left_join(time_means_df, by = "locus_id") %>%
  left_join(time_summary_df, by = "locus_id") %>%
  mutate(
    passed_filter = !is.na(baseMean) & baseMean >= 10,
    filter_status = ifelse(passed_filter, "PASS", "FILTERED_LOW_EXPRESSION"),
    expression_status = ifelse(passed_filter, "Adequate Expression (baseMean >= 10)", "Low/Unexpressed (baseMean < 10)"),
    circadian_rhythmic = passed_filter & !is.na(padj_LRT) & padj_LRT < 0.05,
    rhythmic_category = case_when(
      !passed_filter ~ "Low Expression (baseMean < 10)",
      circadian_rhythmic ~ "Circadian Rhythmic (padj < 0.05)",
      TRUE ~ "Non-Rhythmic (padj >= 0.05)"
    )
  )

write_tsv(thiamine_res, out_tsv)

# 3. 3-panel timecourse figure (scaled to % of peak, 24h wrapping, photoperiod shading)
norm_targets <- norm_counts[targets$locus_id, ]
rownames(norm_targets) <- targets$identifier

replicates_df <- as.data.frame(norm_targets) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "SRR", values_to = "counts") %>%
  left_join(samples %>% select(SRR, time, time_num), by = "SRR")

df_summary <- replicates_df %>%
  group_by(gene, time_num) %>%
  summarise(
    mean_count = mean(counts, na.rm = TRUE),
    sd_count   = sd(counts, na.rm = TRUE),
    sem_count  = sd_count / sqrt(n()),
    .groups = "drop"
  ) %>%
  rename(time = time_num)

df_norm <- df_summary %>%
  group_by(gene) %>%
  mutate(
    peak_val = max(mean_count),
    rel_expr = (mean_count / peak_val) * 100,
    rel_sem  = (sem_count / peak_val) * 100,
    ymin     = pmax(0, rel_expr - rel_sem),
    ymax     = pmin(115, rel_expr + rel_sem)
  ) %>%
  ungroup()

wrap_24h <- df_norm %>% filter(time == 0) %>% mutate(time = 24)
df_plot_norm <- bind_rows(df_norm, wrap_24h) %>% arrange(gene, time)

sunrise_x <- 5.43
sunset_x  <- 19.47

make_panel <- function(gene_order, color_map, shape_map, is_bottom = FALSE, show_labels = FALSE) {
  df_sub <- df_plot_norm %>%
    filter(gene %in% gene_order) %>%
    mutate(gene = factor(gene, levels = gene_order))
  
  p <- ggplot(df_sub, aes(x = time, y = rel_expr, group = gene, color = gene)) +
    annotate("rect", xmin = 0, xmax = sunrise_x, ymin = -Inf, ymax = Inf, fill = "gray92", alpha = 0.6) +
    annotate("rect", xmin = sunset_x, xmax = 24, ymin = -Inf, ymax = Inf, fill = "gray92", alpha = 0.6) +
    geom_vline(xintercept = sunrise_x, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_vline(xintercept = sunset_x,  linetype = "dashed", color = "gray60", linewidth = 0.4)
  
  if (show_labels) {
    p <- p +
      annotate("text", x = 2.7,  y = 108, label = "Dark",  size = 3.2, color = "gray40", fontface = "bold") +
      annotate("text", x = 12.5, y = 108, label = "Light", size = 3.2, color = "gray30", fontface = "bold") +
      annotate("text", x = 21.7, y = 108, label = "Dark",  size = 3.2, color = "gray40", fontface = "bold")
  }
  
  p <- p +
    geom_errorbar(aes(ymin = ymin, ymax = ymax, color = gene), width = 0.5, linewidth = 0.55, alpha = 0.85) +
    geom_line(aes(color = gene), linewidth = 1.05) +
    geom_point(aes(color = gene, shape = gene), size = 2.8, fill = "white", stroke = 1.1) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    scale_x_continuous(breaks = seq(0, 24, 4), labels = c("00:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00"), limits = c(0, 24), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"), expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = if (is_bottom) "Clock Time (Hours)" else NULL, y = "Relative Expression (% of Peak)") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid.major.y = element_line(color = "gray92", linetype = "dotted", linewidth = 0.5),
      axis.title.y = element_text(face = "bold", size = 10),
      axis.title.x = if (is_bottom) element_text(face = "bold", size = 10.5) else element_blank(),
      axis.text.x  = if (is_bottom) element_text(color = "black", size = 9) else element_text(color = "transparent", size = 9),
      axis.ticks.x = if (is_bottom) element_line(color = "black") else element_line(color = "gray80"),
      axis.text.y  = element_text(color = "black", size = 9),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9.5, face = "italic"),
      legend.box.spacing = unit(0.1, "cm"),
      legend.margin = margin(t = 1, b = 2),
      plot.margin = margin(t = 4, r = 12, b = if (is_bottom) 6 else 2, l = 8)
    )
  return(p)
}

panelA <- make_panel(
  gene_order = c("OsCCA1", "OsDof", "OsPCL1"),
  color_map  = c("OsCCA1" = "#D55E00", "OsDof" = "#CC79A7", "OsPCL1" = "#0072B2"),
  shape_map  = c("OsCCA1" = 21, "OsDof" = 24, "OsPCL1" = 22),
  is_bottom  = FALSE,
  show_labels = TRUE
)

panelB <- make_panel(
  gene_order = c("OsDXS1", "OsTHI1", "OsTHIC"),
  color_map  = c("OsDXS1" = "#E69F00", "OsTHI1" = "#CC79A7", "OsTHIC" = "#0072B2"),
  shape_map  = c("OsDXS1" = 21, "OsTHI1" = 22, "OsTHIC" = 24),
  is_bottom  = FALSE,
  show_labels = FALSE
)

panelC <- make_panel(
  gene_order = c("OsGLK1", "OsPIL13"),
  color_map  = c("OsGLK1" = "#009E73", "OsPIL13" = "#D55E00"),
  shape_map  = c("OsGLK1" = 21, "OsPIL13" = 24),
  is_bottom  = TRUE,
  show_labels = FALSE
)


lbl1 <- ggdraw() + draw_label("(I)",   fontfamily = "serif", fontface = "plain", size = 18, x = 0.5, y = 0.55)
lbl2 <- ggdraw() + draw_label("(II)",  fontfamily = "serif", fontface = "plain", size = 18, x = 0.5, y = 0.55)
lbl3 <- ggdraw() + draw_label("(III)", fontfamily = "serif", fontface = "plain", size = 18, x = 0.5, y = 0.59)

row1 <- plot_grid(lbl1, panelA, ncol = 2, rel_widths = c(0.10, 0.90))
row2 <- plot_grid(lbl2, panelB, ncol = 2, rel_widths = c(0.10, 0.90))
row3 <- plot_grid(lbl3, panelC, ncol = 2, rel_widths = c(0.10, 0.90))

p_3panel <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1, 1, 1))

for (d in c(pub_dir, out_dir, thesis_fig_dir)) {
  ggsave(file.path(d, "main_fig_timecourse_3panels.pdf"), plot = p_3panel, width = 5.2, height = 9.5, device = cairo_pdf)
  ggsave(file.path(d, "main_fig_timecourse_3panels.png"), plot = p_3panel, width = 5.2, height = 9.5, dpi = 300)
}

# 4. pairwise Pearson correlation heatmap on VST vectors (complete linkage, Euclidean distance)
expressed_targets <- thiamine_res %>% filter(passed_filter == TRUE)

vst_expressed <- vst_mat[expressed_targets$locus_id, ]
rownames(vst_expressed) <- expressed_targets$identifier

cor_matrix <- cor(t(vst_expressed), method = "pearson")
cor_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
breaks_seq <- seq(-1, 1, length.out = 101)

pdf(file.path(pub_dir, "annex_correlation_matrix.pdf"), width = 6.8, height = 9.5)
pheatmap(
  cor_matrix,
  color = cor_palette,
  breaks = breaks_seq,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 9.5,
  fontsize_row = 9.0,
  fontsize_col = 9.0,
  treeheight_row = 35,
  treeheight_col = 35,
  border_color = "gray85",
  main = NA
)
invisible(dev.off())

png(file.path(pub_dir, "annex_correlation_matrix.png"), width = 6.8, height = 9.5, units = "in", res = 300)
pheatmap(
  cor_matrix,
  color = cor_palette,
  breaks = breaks_seq,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 9.5,
  fontsize_row = 9.0,
  fontsize_col = 9.0,
  treeheight_row = 35,
  treeheight_col = 35,
  border_color = "gray85",
  main = NA
)
invisible(dev.off())
