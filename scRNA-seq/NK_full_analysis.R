# =============================================================================
# NK_full_analysis.R — NK cell main analysis: TIGIT expression, metabolic
#                      pathways, cytotoxicity/exhaustion modules
# TB-TIGIT-NK manuscript
#
# Input:
#   ./data/mtb_final.rds          (Seurat, ~358K PBMC, 120 donors)
#   ./data/Lung_HC_TB_annotated_0.5.rds  (Seurat, for Part 8 module scores)
#   ./data/114条代谢通路基因.xlsx          (Excel, for Parts 2/5/6)
#
# Output panels:
#   1. TIGIT expression distribution (histogram + paired dot plot)
#   2. 114 metabolic pathway heatmap (Kruskal-Wallis ranked)     [if(FALSE)]
#   3. TIGIT Low/High grouping histogram
#   4. DotPlot: effector genes in TIGIT Low vs High
#   5. 24h TIGIT Low/High — 114 pathway heatmap                  [if(FALSE)]
#   6. 24h TIGIT Low/High — key metabolomics modules heatmap     [if(FALSE)]
#   7. TIGIT Low vs High boxplot
#   8. Cytotoxicity/Exhaustion donor-level module scores
# =============================================================================

# ===== SETUP =====
rm(list = ls())
cat("=== NK Full Analysis ===\n")

library(Seurat)
library(ggplot2)
library(pheatmap)
library(readxl)
library(dplyr)
library(tidyr)
library(ggpubr)
library(patchwork)
library(gghalves)

source("utils/scRNA-seq.R")

# Output and data directories
OUT_DIR <- "./output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
DATA_DIR <- "./data"

set.seed(42)

# ===== LOAD DATA =====
cat("\n===== Loading data =====\n")
obj <- readRDS(file.path(DATA_DIR, "mtb_final.rds"))
nk <- obj[, obj$cell_type_lowerres == "NK"]
nk$timepoint <- factor(nk$timepoint, levels = c("UT", "3hMTB", "24hMTB"))
cat(sprintf("NK cells: %d\n", ncol(nk)))
print(table(nk$timepoint))

tigit_expr <- GetAssayData(nk, layer = "data")["TIGIT", ]
col_tp <- c("UT" = "#E69F00", "3hMTB" = "#D55E00", "24hMTB" = "#CC79A7")

# =============================================================================
# PART 1: TIGIT Expression Distribution
# =============================================================================
cat("\n===== PART 1: TIGIT Expression =====\n")

tig_summary <- data.frame()
for (tp in c("UT", "3hMTB", "24hMTB")) {
  cells_tp <- colnames(nk)[nk$timepoint == tp]
  e <- tigit_expr[cells_tp]
  tig_summary <- rbind(tig_summary, data.frame(
    timepoint = tp, total = length(e),
    pos = sum(e > 0), pct = mean(e > 0) * 100,
    median_pos = median(e[e > 0])
  ))
}
print(tig_summary)

# Histogram (TIGIT+ only)
all_pos <- tigit_expr[tigit_expr > 0]
df_tig <- data.frame(
  expr = all_pos,
  timepoint = nk[, names(all_pos)]$timepoint
)

p1 <- ggplot(df_tig, aes(x = expr, fill = timepoint)) +
  geom_histogram(bins = 60, alpha = 0.7, color = "white", linewidth = 0.15,
                 position = "identity") +
  scale_fill_manual(values = col_tp) +
  labs(title = "TIGIT Expression in NK cells (TIGIT+ only)",
       x = "TIGIT Expression (log-normalized)", y = "Cell Count",
       fill = "Timepoint") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "top")

# Donor-level paired dot plot
pd <- position_jitter(width = 0.1, seed = 42)
my_comp <- list(c("UT", "3hMTB"), c("UT", "24hMTB"))

nk_avg <- FetchData(nk, vars = c("TIGIT", "timepoint", "assignment")) %>%
  filter(assignment != "") %>%
  group_by(assignment, timepoint) %>%
  summarise(TIGIT = mean(TIGIT, na.rm = TRUE), .groups = "drop") %>%
  group_by(assignment) %>%
  filter(n() == 3) %>%
  ungroup() %>%
  mutate(timepoint = factor(timepoint, levels = c("UT", "3hMTB", "24hMTB")))

p2 <- ggplot(nk_avg, aes(x = timepoint, y = TIGIT)) +
  geom_boxplot(aes(fill = timepoint), outlier.shape = NA, alpha = 0.6, width = 0.5) +
  geom_line(aes(group = assignment), color = "gray60", alpha = 0.2,
            linetype = "twodash", position = pd) +
  geom_point(size = 1.5, shape = 21, fill = "white", color = "black",
             stroke = 0.8, position = pd) +
  scale_fill_manual(values = col_tp) +
  stat_compare_means(comparisons = my_comp, paired = TRUE,
                     method = "wilcox.test", label = "p.signif",
                     tip.length = 0.01, size = 5) +
  labs(title = "TIGIT Expression in NK cells (Paired Donors)",
       y = "Average Expression", x = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "none")

p_combined <- p1 | p2
ggsave(file.path(OUT_DIR, "NK_TIGIT_expression.pdf"), p_combined, width = 8, height = 5, dpi = 300)
cat("  Saved: NK_TIGIT_expression.pdf\n")

# =============================================================================
# PART 2: 114 Metabolic Pathway Heatmap (Kruskal-Wallis ranked)
# =============================================================================
cat("\n===== PART 2: 114 Metabolic Pathway Heatmap (SKIP) =====\n")

if (FALSE) {
  cat("  >> Running 114-pathway heatmap...\n")

  df <- read_excel(file.path(DATA_DIR, "114\u6761\u4ee3\u8c22\u901a\u8def\u57fa\u56e0.xlsx"), sheet = "Sheet1")
  pathway_list <- list()
  for (col in colnames(df)) {
    genes <- df[[col]][!is.na(df[[col]])]
    genes <- trimws(genes); genes <- genes[genes != ""]
    if (length(genes) >= 3) {
      g <- intersect(genes, rownames(nk))
      if (length(g) >= 3) pathway_list[[col]] <- g
    }
  }
  cat(sprintf("  Valid pathways: %d\n", length(pathway_list)))

  expr_mat <- GetAssayData(nk, layer = "data")
  score_mat <- sapply(pathway_list, function(genes) {
    colMeans(expr_mat[genes, , drop = FALSE])
  })
  rownames(score_mat) <- colnames(nk)

  meta_df <- data.frame(timepoint = nk$timepoint, row.names = colnames(nk))
  summ <- score_mat %>%
    as.data.frame() %>%
    cbind(meta_df) %>%
    pivot_longer(-timepoint, names_to = "pathway", values_to = "score") %>%
    group_by(pathway, timepoint) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = timepoint, values_from = mean_score) %>%
    as.data.frame()
  rownames(summ) <- summ$pathway; summ$pathway <- NULL

  heat_mat <- t(scale(t(as.matrix(summ))))

  kw_p <- apply(score_mat, 2, function(x) {
    kruskal.test(x ~ nk$timepoint)$p.value
  })
  kw_df <- data.frame(pathway = names(kw_p), p_value = kw_p)
  kw_df <- kw_df[order(kw_df$p_value), ]

  heat_plot <- heat_mat[kw_df$pathway, ]
  max_abs <- max(abs(heat_plot))

  timepoint_colors <- c("UT" = "#E69F00", "3hMTB" = "#D55E00", "24hMTB" = "#CC79A7")
  annotation_col <- data.frame(
    Timepoint = factor(colnames(heat_plot), levels = c("UT", "3hMTB", "24hMTB")),
    row.names = colnames(heat_plot)
  )

  pdf(file.path(OUT_DIR, "NK_114pathways_heatmap.pdf"), width = 8, height = 22)
  pheatmap(heat_plot,
           color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
           breaks = seq(-max_abs, max_abs, length.out = 101),
           cluster_rows = TRUE, cluster_cols = FALSE,
           annotation_col = annotation_col, cellwidth = 27, cellheight = 7.5,
           annotation_colors = list(Timepoint = timepoint_colors),
           fontsize_row = 7, fontsize_col = 12,
           main = "Metabolic Pathways in NK cells (All 114)",
           border_color = NA)
  dev.off()
  cat("  Saved: NK_114pathways_heatmap.pdf\n")
}

# =============================================================================
# PART 3: TIGIT Low/High Grouping Histogram
# =============================================================================
cat("\n===== PART 3: TIGIT Low/High Grouping =====\n")

n <- 800
all_low <- c(); all_high <- c()

hist_list <- list()
for (tp in c("UT", "3hMTB", "24hMTB")) {
  cells_tp <- colnames(nk)[nk$timepoint == tp]
  tig_tp <- tigit_expr[cells_tp]
  tig_pos <- tig_tp[tig_tp > 0]
  tig_sorted <- sort(tig_pos)

  low_cells <- names(tig_sorted[1:n])
  high_cells <- names(rev(tig_sorted)[1:n])
  low_max <- tig_sorted[n]
  high_min <- rev(tig_sorted)[n]

  all_low <- c(all_low, low_cells)
  all_high <- c(all_high, high_cells)

  cat(sprintf("  %s: Low=%.4f~%.4f, High=%.4f~%.4f\n",
      tp, min(tig_sorted[1:n]), low_max, high_min, max(rev(tig_sorted)[1:n])))

  df_h <- data.frame(expr = tig_pos, timepoint = tp)
  df_h$group <- case_when(
    df_h$expr <= low_max ~ "Low (selected)",
    df_h$expr >= high_min ~ "High (selected)",
    TRUE ~ "Middle (excluded)"
  )
  df_h$group <- factor(df_h$group, levels = c("Low (selected)", "Middle (excluded)", "High (selected)"))
  hist_list[[tp]] <- df_h
}

df_hist <- do.call(rbind, hist_list)
df_hist$timepoint <- factor(df_hist$timepoint, levels = c("UT", "3hMTB", "24hMTB"))

p_hist <- ggplot(df_hist, aes(x = expr, fill = group)) +
  geom_histogram(bins = 50, alpha = 0.85, color = "white", linewidth = 0.12) +
  scale_fill_manual(values = c("Low (selected)" = "#E8A838",
                                "Middle (excluded)" = "grey80",
                                "High (selected)" = "#C43A3A"),
                    guide = guide_legend(direction = "horizontal")) +
  facet_wrap(~ timepoint, ncol = 3, scales = "free_y") +
  labs(title = sprintf("TIGIT Low vs High in NK cells (n=%d each)", n),
       x = "TIGIT Expression", y = "Count") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.position = "top",
        legend.margin = margin(b = -5),
        plot.margin = margin(5, 5, 5, 5))

ggsave(file.path(OUT_DIR, "NK_TIGIT_group_hist.pdf"), p_hist, width = 9, height = 3.5, dpi = 300)
cat("  Saved: NK_TIGIT_group_hist.pdf\n")

# =============================================================================
# PART 4: DotPlot — Effector Genes in TIGIT Low vs High
# =============================================================================
cat("\n===== PART 4: DotPlot =====\n")

cells_use <- c(all_low, all_high)
nk_sub <- nk[, cells_use]
nk_sub$TIGIT_group <- ifelse(colnames(nk_sub) %in% all_high, "High", "Low")
nk_sub$group_label <- paste0(nk_sub$timepoint, "_", nk_sub$TIGIT_group)
nk_sub$group_label <- factor(nk_sub$group_label,
  levels = c("UT_Low", "UT_High", "3hMTB_Low", "3hMTB_High", "24hMTB_Low", "24hMTB_High"))

cat("  Group counts:\n")
print(table(nk_sub$group_label))

genes_dot <- c("TNF", "TNFSF10", "PRF1", "CCL5", "GZMA",
               "GZMB", "IFNG", "IL1B", "IL6", "LAMP1")

genes_avail <- intersect(genes_dot, rownames(nk_sub))
cat(sprintf("  DotPlot genes: %d/%d\n", length(genes_avail), length(genes_dot)))
if (length(genes_avail) < length(genes_dot)) {
  cat("  Missing:", setdiff(genes_dot, genes_avail), "\n")
}

Idents(nk_sub) <- "group_label"

p_dot <- DotPlot(nk_sub, features = genes_avail,
                 cols = c("lightgrey", "#C43A3A"), dot.scale = 10) +
  labs(title = "TIGIT Low vs High in NK cells (n=800 per group per timepoint)",
       x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10))

ggsave(file.path(OUT_DIR, "NK_TIGIT_dotplot.pdf"), p_dot, width = 8, height = 4, dpi = 300)
cat("  Saved: NK_TIGIT_dotplot.pdf\n")

# =============================================================================
# PART 5: 24h TIGIT Low/High — 114 Pathway Heatmap
# =============================================================================
cat("\n===== PART 5: 24h TIGIT Low/High 114-pathway Heatmap (SKIP) =====\n")

if (FALSE) {
  cat("  >> Running 24h TIGIT Low/High 114-pathway heatmap...\n")

  nk_24h <- nk[, nk$timepoint == "24hMTB"]
  tigit_expr_24h <- GetAssayData(nk_24h, layer = "data")["TIGIT", ]
  tigit_pos_24h <- tigit_expr_24h[tigit_expr_24h > 0]
  tigit_sorted_24h <- sort(tigit_pos_24h)

  n_24 <- 800
  low_cells_24 <- names(tigit_sorted_24h[1:n_24])
  high_cells_24 <- names(rev(tigit_sorted_24h)[1:n_24])
  nk_sub_24 <- nk_24h[, c(low_cells_24, high_cells_24)]
  nk_sub_24$TIGIT_group <- factor(ifelse(colnames(nk_sub_24) %in% high_cells_24, "High", "Low"),
                                   levels = c("Low", "High"))

  # Read and filter pathways
  df_p5 <- read_excel(file.path(DATA_DIR, "114\u6761\u4ee3\u8c22\u901a\u8def\u57fa\u56e0.xlsx"), sheet = "Sheet1")
  pathway_list_p5 <- list()
  for (col in colnames(df_p5)) {
    genes <- df_p5[[col]][!is.na(df_p5[[col]])]
    genes <- trimws(genes); genes <- genes[genes != ""]
    if (length(genes) >= 3) {
      g <- intersect(genes, rownames(nk_sub_24))
      if (length(g) >= 3) pathway_list_p5[[col]] <- g
    }
  }
  cat(sprintf("  Valid pathways: %d\n", length(pathway_list_p5)))

  # Pathway scores + z-score
  expr_mat_p5 <- GetAssayData(nk_sub_24, layer = "data")
  score_mat_p5 <- sapply(pathway_list_p5, function(gs) colMeans(expr_mat_p5[gs, , drop = FALSE]))
  rownames(score_mat_p5) <- colnames(nk_sub_24)
  score_mat_z_p5 <- scale(score_mat_p5)

  # Mean by TIGIT group
  meta_p5 <- data.frame(TIGIT_group = nk_sub_24$TIGIT_group, row.names = colnames(nk_sub_24))
  summ_p5 <- score_mat_z_p5 %>%
    as.data.frame() %>%
    cbind(meta_p5) %>%
    pivot_longer(-TIGIT_group, names_to = "pathway", values_to = "score") %>%
    group_by(pathway, TIGIT_group) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = TIGIT_group, values_from = mean_score) %>%
    as.data.frame()
  rownames(summ_p5) <- summ_p5$pathway; summ_p5$pathway <- NULL

  # Wilcoxon ranking
  wilcox_p_p5 <- apply(score_mat_z_p5, 2, function(x) {
    wilcox.test(x ~ nk_sub_24$TIGIT_group, exact = FALSE)$p.value
  })
  kw_df_p5 <- data.frame(pathway = names(wilcox_p_p5), p_value = wilcox_p_p5)
  kw_df_p5 <- kw_df_p5[order(kw_df_p5$p_value), ]

  heat_plot_p5 <- as.matrix(summ_p5[kw_df_p5$pathway, ])
  max_abs_p5 <- max(abs(heat_plot_p5))

  annotation_col_p5 <- data.frame(
    Group = factor(colnames(heat_plot_p5), levels = c("Low", "High")),
    row.names = colnames(heat_plot_p5)
  )
  annotation_colors_p5 <- list(Group = c("Low" = "#E8A838", "High" = "#C43A3A"))

  pdf(file.path(OUT_DIR, "NK_24h_TIGIT_metabolic_heatmap.pdf"), width = 5, height = 22)
  pheatmap(heat_plot_p5,
           color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
           breaks = seq(-max_abs_p5, max_abs_p5, length.out = 101),
           cluster_rows = TRUE, cluster_cols = FALSE,
           annotation_col = annotation_col_p5, cellwidth = 24, cellheight = 7.5,
           annotation_colors = annotation_colors_p5,
           fontsize_row = 7, fontsize_col = 12,
           main = "24h NK: TIGIT Low vs High (114 Metabolic Pathways)",
           border_color = NA)
  dev.off()
  cat("  Saved: NK_24h_TIGIT_metabolic_heatmap.pdf\n")
}

# =============================================================================
# PART 6: 24h TIGIT Low/High — Key Metabolomics Modules Heatmap
# =============================================================================
cat("\n===== PART 6: 24h Key Metabolomics Modules Heatmap (SKIP) =====\n")

if (FALSE) {
  cat("  >> Running 24h key metabolomics modules heatmap...\n")

  nk_24h_p6 <- nk[, nk$timepoint == "24hMTB"]
  tigit_expr_p6 <- GetAssayData(nk_24h_p6, layer = "data")["TIGIT", ]
  tigit_pos_p6 <- tigit_expr_p6[tigit_expr_p6 > 0]
  tigit_sorted_p6 <- sort(tigit_pos_p6)

  n_p6 <- 800
  low_cells_p6 <- names(tigit_sorted_p6[1:n_p6])
  high_cells_p6 <- names(rev(tigit_sorted_p6)[1:n_p6])
  cells_use_p6 <- c(low_cells_p6, high_cells_p6)

  nk_sub_p6 <- nk_24h_p6[, cells_use_p6]
  nk_sub_p6$TIGIT_group <- factor(ifelse(colnames(nk_sub_p6) %in% high_cells_p6, "High", "Low"),
                                   levels = c("Low", "High"))

  # Key pathways
  df_p6 <- read_excel(file.path(DATA_DIR, "114\u6761\u4ee3\u8c22\u901a\u8def\u57fa\u56e0.xlsx"), sheet = "Sheet1")
  pathway_map <- list(
    Pyrimidine   = "Pyrimidine Metabolism",
    CoA          = "Pantothenate and CoA Biosynthesis",
    Polyamine    = "Polyamine Biosynthesis",
    OXPHOS       = grep("Oxidative", colnames(df_p6), value = TRUE)[1],
    TCA          = "Citric Acid Cycle",
    Glycolysis   = "Glycolysis",
    Tryptophan   = "Tryptophan Metabolism",
    Kynurenine   = "Kynurenine Metabolism",
    Beta_Alanine = "Beta-Alanine Metabolism",
    Glutamate    = "Alanine, Aspartate and Glutamate Metabolism"
  )
  pathway_list_p6 <- list()
  for (nm in names(pathway_map)) {
    cn <- pathway_map[[nm]]
    if (cn %in% colnames(df_p6)) {
      g <- intersect(trimws(na.omit(df_p6[[cn]])), rownames(nk_sub_p6))
      if (length(g) >= 3) pathway_list_p6[[nm]] <- g
    }
  }

  # colMeans + z-score
  expr_p6 <- GetAssayData(nk_sub_p6, layer = "data")
  score_mat_p6 <- sapply(pathway_list_p6, function(g) colMeans(expr_p6[g, , drop = FALSE]))
  score_z_p6 <- scale(score_mat_p6)

  # Statistics
  meta_p6 <- data.frame(TIGIT_group = nk_sub_p6$TIGIT_group, score_z_p6)
  stats_p6 <- meta_p6 %>%
    group_by(TIGIT_group) %>%
    summarise(across(everything(), mean)) %>%
    as.data.frame()
  rownames(stats_p6) <- stats_p6$TIGIT_group; stats_p6$TIGIT_group <- NULL

  wilcox_list_p6 <- list()
  for (nm in names(pathway_list_p6)) {
    w <- wilcox.test(meta_p6[[nm]] ~ meta_p6$TIGIT_group, exact = FALSE)
    wilcox_list_p6[[nm]] <- data.frame(
      pathway = nm,
      Low = stats_p6["Low", nm], High = stats_p6["High", nm],
      diff = stats_p6["High", nm] - stats_p6["Low", nm],
      p_value = w$p.value, stringsAsFactors = FALSE
    )
  }
  result_p6 <- bind_rows(wilcox_list_p6)
  result_p6$p_adj <- p.adjust(result_p6$p_value, "BH")
  result_p6$sig <- ifelse(result_p6$p_adj < 0.001, "***",
                   ifelse(result_p6$p_adj < 0.01, "**",
                   ifelse(result_p6$p_adj < 0.05, "*", "ns")))
  result_p6 <- result_p6[order(abs(result_p6$diff), decreasing = TRUE), ]
  cat("  === Results (colMeans + z-score) ===\n")
  print(result_p6)

  # Vertical heatmap
  heat_p6 <- as.matrix(result_p6[, c("Low", "High")])
  rownames(heat_p6) <- result_p6$pathway
  max_abs_p6 <- max(abs(heat_p6))

  annotation_col_p6 <- data.frame(Group = c("Low", "High"), row.names = c("Low", "High"))

  pdf(file.path(OUT_DIR, "NK_24h_metabolomics_modules_v2.pdf"), width = 4, height = 5)
  pheatmap(heat_p6,
           color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
           breaks = seq(-max_abs_p6, max_abs_p6, length.out = 101),
           cluster_rows = TRUE, cluster_cols = FALSE,
           annotation_col = annotation_col_p6, cellwidth = 25, cellheight = 20,
           annotation_colors = list(Group = c("Low" = "#E8A838", "High" = "#C43A3A")),
           display_numbers = FALSE, number_format = "%.2f",
           fontsize_row = 9, fontsize_col = 11,
           main = "24h NK: TIGIT Low vs High\nMetabolic Pathways (z-score, 2-group)",
           border_color = NA)
  dev.off()
  cat("  Saved: NK_24h_metabolomics_modules_v2.pdf\n")

  # Horizontal heatmap (Low/High as rows, pathways as columns)
  heat_t_p6 <- t(heat_p6)
  max_abs_t_p6 <- max(abs(heat_t_p6))

  annotation_row_p6 <- data.frame(Group = c("Low", "High"), row.names = c("Low", "High"))

  pdf(file.path(OUT_DIR, "NK_24h_metabolomics_modules_v2_horizontal.pdf"), width = 7, height = 2.5)
  pheatmap(heat_t_p6,
           color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
           breaks = seq(-max_abs_t_p6, max_abs_t_p6, length.out = 101),
           cluster_rows = FALSE, cluster_cols = TRUE,
           annotation_row = annotation_row_p6, cellwidth = 20, cellheight = 22,
           annotation_colors = list(Group = c("Low" = "#E8A838", "High" = "#C43A3A")),
           display_numbers = FALSE, number_format = "%.2f",
           fontsize_row = 11, fontsize_col = 10,
           main = "24h NK: TIGIT Low vs High",
           border_color = NA, angle_col = 45)
  dev.off()
  cat("  Saved: NK_24h_metabolomics_modules_v2_horizontal.pdf\n")
}

# =============================================================================
# PART 7: TIGIT Low vs High Boxplot
# =============================================================================
cat("\n===== PART 7: TIGIT Low vs High Boxplot =====\n")

# Re-use all_low / all_high from Part 3
cells_use_p7 <- c(all_low, all_high)

df_plot_p7 <- data.frame(
  cell = cells_use_p7,
  TIGIT_expr = tigit_expr[cells_use_p7],
  group = rep(c("TIGIT_Low", "TIGIT_High"), each = n)
)
df_plot_p7$group <- factor(df_plot_p7$group, levels = c("TIGIT_Low", "TIGIT_High"))

p_box <- ggplot(df_plot_p7, aes(x = group, y = TIGIT_expr, fill = group, color = group)) +
  geom_jitter(width = 0.12, size = 0.4, alpha = 0.35) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black", width = 0.4) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 8) +
  scale_fill_manual(values = c("TIGIT_Low" = "#E8A838", "TIGIT_High" = "#C43A3A")) +
  scale_color_manual(values = c("TIGIT_Low" = "#E8A838", "TIGIT_High" = "#C43A3A")) +
  labs(title = "TIGIT Expression: Low vs High",
       x = NULL, y = "TIGIT Expression") +
  theme_classic() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

ggsave(file.path(OUT_DIR, "TIGIT_low_high_boxplot.pdf"), p_box, width = 6, height = 7, dpi = 300)
cat("  Saved: TIGIT_low_high_boxplot.pdf\n")

# =============================================================================
# PART 8: Cytotoxicity / Exhaustion Donor-Level Module Scores
# =============================================================================
cat("\n===== PART 8: Cytotoxicity / Exhaustion Module Scores =====\n")

obj_p8 <- readRDS(file.path(DATA_DIR, "Lung_HC_TB_annotated_0.5.rds"))
nk_p8 <- subset(obj_p8, idents = "NK")
nk_p8$group <- factor(nk_p8$group, levels = c("HC", "FDR-LOW", "FDR-HIGH"))

module_genes <- list(
  Cytotoxicity = c("GZMA","GZMB","PRF1","GNLY","FASLG","NKG7","KLRD1","LAMP1"),
  Exhaustion   = c("PDCD1","CTLA4","HAVCR2","LAG3","TOX","KLRC1")
)

for (nm in names(module_genes)) {
  g <- intersect(module_genes[[nm]], rownames(nk_p8))
  module_genes[[nm]] <- g
  cat(sprintf("  %s: %d genes\n", nm, length(g)))
}

for (nm in names(module_genes)) {
  nk_p8 <- AddModuleScore(nk_p8, features = list(module_genes[[nm]]), name = "tmp_", seed = 42)
  colnames(nk_p8@meta.data)[ncol(nk_p8@meta.data)] <- nm
}

donor_p8 <- nk_p8@meta.data %>%
  group_by(orig.ident, group) %>%
  summarise(Cyto = mean(Cytotoxicity), Exh = mean(Exhaustion), .groups = "drop")

donor_long_p8 <- donor_p8 %>%
  pivot_longer(cols = c("Cyto", "Exh"), names_to = "Module", values_to = "Score")

make_plot_p8 <- function(data, score_col, title, label_y) {
  ggplot(data, aes(x = group, y = .data[[score_col]], fill = group)) +
    gghalves::geom_half_violin(side = "r", alpha = 0.5, scale = "width", trim = FALSE) +
    geom_jitter(width = 0.06, size = 3, shape = 21, color = "brown", stroke = 0.4, alpha = 0.8) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.7) +
    stat_compare_means(
      comparisons = list(c("HC","FDR-LOW"), c("HC","FDR-HIGH"), c("FDR-LOW","FDR-HIGH")),
      method = "wilcox.test", label = "p.signif", size = 6,
      label.y = label_y) +
    scale_fill_manual(values = c("HC"="#374e55","FDR-LOW"="#7eb36e","FDR-HIGH"="#c44e52")) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          legend.position = "none") +
    labs(title = title, x = "", y = "Module Score")
}

pp1 <- make_plot_p8(donor_long_p8[donor_long_p8$Module == "Cyto", ], "Score", "Cytotoxicity", c(1.2, 1.3, 1.4))
pp2 <- make_plot_p8(donor_long_p8[donor_long_p8$Module == "Exh", ], "Score", "Exhaustion", c(0.07, 0.09, 0.11))

pp <- pp1 | pp2
ggsave(file.path(OUT_DIR, "NK_FDR_donor_2_module_scores.pdf"), plot = pp, width = 10, height = 5, dpi = 300)
cat("  Saved: NK_FDR_donor_2_module_scores.pdf\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n========================================\n")
cat("All outputs saved to:", normalizePath(OUT_DIR), "\n")
cat("========================================\n")
cat("  Panel 1: NK_TIGIT_expression.pdf\n")
cat("  Panel 2: NK_114pathways_heatmap.pdf          [if(FALSE)]\n")
cat("  Panel 3: NK_TIGIT_group_hist.pdf\n")
cat("  Panel 4: NK_TIGIT_dotplot.pdf\n")
cat("  Panel 5: NK_24h_TIGIT_metabolic_heatmap.pdf  [if(FALSE)]\n")
cat("  Panel 6: NK_24h_metabolomics_modules_v2.pdf  [if(FALSE)]\n")
cat("            NK_24h_metabolomics_modules_v2_horizontal.pdf  [if(FALSE)]\n")
cat("  Panel 7: TIGIT_low_high_boxplot.pdf\n")
cat("  Panel 8: NK_FDR_donor_2_module_scores.pdf\n")
cat("========================================\n")
cat("=== NK Full Analysis Complete ===\n")
