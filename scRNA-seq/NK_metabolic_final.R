# ============================================================
# SETUP
# ============================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(pheatmap)
library(patchwork)

source("utils/scRNA-seq.R")

set.seed(42)

DATA_DIR <- "./data"
OUT_DIR  <- "./output"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== NK_metabolic_final.R ===\n\n")
cat(sprintf("Data dir : %s\n", DATA_DIR))
cat(sprintf("Output dir: %s\n\n", OUT_DIR))

# ============================================================
# LOAD DATA & SUBSET NK
# ============================================================
cat("Loading data ...\n")
obj <- readRDS(file.path(DATA_DIR, "mtb_final.rds"))
nk  <- obj[, obj$cell_type_lowerres == "NK"]
nk$timepoint <- factor(nk$timepoint, levels = c("UT", "3hMTB", "24hMTB"))
cat(sprintf("NK cells: %d\n", ncol(nk)))

# ============================================================
# COLOR & THEME
# ============================================================
col_tp <- c("UT" = "#E69F00", "3hMTB" = "#D55E00", "24hMTB" = "#CC79A7")

theme_my <- theme_classic() +
  theme(
    axis.text.x    = element_text(color = "black", size = 12, face = "bold"),
    axis.text.y    = element_text(color = "black", size = 11, face = "bold"),
    axis.title     = element_text(face = "bold", size = 13),
    plot.title     = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "none"
  )

# ============================================================
# GENE LISTS
# ============================================================
gly_keep <- c("SLC2A1", "HK2", "PFKFB3", "PKM", "ENO1", "LDHA")
lac_keep <- c("LDHA", "SLC16A1", "SLC16A3", "HIF1A")
oxp_keep <- c("COX5A", "NDUFA1", "NDUFS1", "NDUFS3", "NDUFV1", "TFAM")
all_keep <- unique(c(gly_keep, lac_keep, oxp_keep))

# ============================================================
# DONOR-AVERAGE EXPRESSION
# ============================================================
cat("Computing donor averages ...\n")
meta_cols <- c("assignment", "timepoint")
expr_df <- FetchData(nk, vars = c(meta_cols, all_keep), layer = "data")
donor_avg <- expr_df %>%
  filter(assignment != "") %>%
  group_by(assignment, timepoint) %>%
  summarise(across(all_of(all_keep), mean), .groups = "drop") %>%
  group_by(assignment) %>%
  filter(n() == 3) %>%
  ungroup()

donor_long <- donor_avg %>%
  pivot_longer(
    cols      = all_of(all_keep),
    names_to  = "gene",
    values_to = "expr"
  ) %>%
  mutate(timepoint = factor(timepoint, levels = c("UT", "3hMTB", "24hMTB")))

pd       <- position_jitter(width = 0.1, seed = 42)
my_comp  <- list(c("UT", "3hMTB"), c("UT", "24hMTB"))
cat(sprintf("Donors with complete timepoints: %d\n", n_distinct(donor_avg$assignment)))

# ============================================================
# SECTION 1: INDEPENDENT BOXPLOTS (PER GENE)
# ============================================================
cat("\n=== 1. Boxplots ===\n")

build_gene_plot <- function(df, g) {
  sub <- df %>% filter(gene == g)
  ggplot(sub, aes(x = timepoint, y = expr)) +
    geom_boxplot(aes(fill = timepoint), outlier.shape = NA, alpha = 0.6, width = 0.5) +
    geom_line(aes(group = assignment), color = "gray60", alpha = 0.2,
              linetype = "twodash", position = pd) +
    geom_point(size = 1.5, shape = 21, fill = "white", color = "black",
               stroke = 0.8, position = pd) +
    scale_fill_manual(values = col_tp) +
    stat_compare_means(comparisons = my_comp, paired = TRUE, method = "wilcox.test",
                       label = "p.signif", tip.length = 0.01, size = 5) +
    labs(title = g, y = "Average Expression", x = "") +
    theme_my
}

combine_module <- function(df, genes, title, fname, ncol = 4) {
  plots <- lapply(genes, build_gene_plot, df = df)
  nrow_plot <- ceiling(length(plots) / ncol)
  w <- ncol * 3.5 + 1
  h <- nrow_plot * 4 + 0.8
  p <- wrap_plots(plots, ncol = ncol) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
    )
  ggsave(file.path(OUT_DIR, fname), p, width = w, height = h, dpi = 300)
  cat(sprintf("  %s (%d genes, %dx%d, %.1f x %.1f)\n",
              fname, length(genes), ncol, nrow_plot, w, h))
}

combine_module(donor_long, gly_keep,
               "Glycolysis Genes in NK cells",
               "NK_Glycolysis_PairedBox.pdf", 4)
combine_module(donor_long, lac_keep,
               "Lactate Metabolism Genes in NK cells",
               "NK_Lactate_PairedBox.pdf", 4)
combine_module(donor_long, oxp_keep,
               "OXPHOS/Mitochondrial Genes in NK cells",
               "NK_OXPHOS_PairedBox.pdf", 3)

# ============================================================
# SECTION 2: TIMECOURSE HEATMAP
# ============================================================
cat("\n=== 2. Timecourse heatmap ===\n")

tp_means <- donor_avg %>%
  group_by(timepoint) %>%
  summarise(across(all_of(all_keep), mean), .groups = "drop") %>%
  mutate(timepoint = factor(timepoint, levels = c("UT", "3hMTB", "24hMTB"))) %>%
  arrange(timepoint) %>%
  as.data.frame()
rownames(tp_means) <- tp_means$timepoint
tp_means$timepoint <- NULL
tp_means <- tp_means[c("UT", "3hMTB", "24hMTB"), ]

heat1 <- t(scale(as.matrix(tp_means)))
heat1[heat1 > 2] <- 2
heat1[heat1 < -2] <- -2
heat1 <- heat1[, c("UT", "3hMTB", "24hMTB"), drop = FALSE]
max_abs <- max(abs(heat1), na.rm = TRUE)

ann_col <- data.frame(
  Timepoint = factor(c("UT", "3hMTB", "24hMTB"),
                     levels = c("UT", "3hMTB", "24hMTB")),
  row.names = c("UT", "3hMTB", "24hMTB")
)

pdf(file.path(OUT_DIR, "NK_Metabolic_Timecourse_Heatmap.pdf"),
    width = 6, height = 8)
pheatmap(
  heat1,
  color             = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  breaks            = seq(-max_abs, max_abs, length.out = 101),
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  annotation_col    = ann_col,
  cellwidth         = 26,
  cellheight        = 14,
  annotation_colors = list(Timepoint = col_tp),
  fontsize_row      = 11,
  fontsize_col      = 12,
  border_color      = NA,
  main              = "NK Metabolic Genes across Timepoints"
)
dev.off()
cat("  Done\n")

# ============================================================
# SECTION 3: TIGIT HIGH vs LOW (3h + 24h)
# ============================================================
cat("\n=== 3. TIGIT High vs Low ===\n")

process_tigit <- function(tp, n_each = 800) {
  nk_tp <- nk[, nk$timepoint == tp]
  tigit_expr <- GetAssayData(nk_tp, layer = "data")["TIGIT", ]
  tigit_pos <- tigit_expr[tigit_expr > 0]
  tigit_sorted <- sort(tigit_pos)

  if (length(tigit_sorted) < n_each * 2) {
    n_each <- floor(length(tigit_sorted) / 2)
  }
  low_cells  <- names(tigit_sorted[1:n_each])
  high_cells <- names(rev(tigit_sorted)[1:n_each])

  nk_sub <- nk_tp[, c(low_cells, high_cells)]
  nk_sub$TIGIT_group <- factor(
    ifelse(colnames(nk_sub) %in% high_cells, "High", "Low"),
    levels = c("Low", "High")
  )
  list(obj = nk_sub, n = n_each)
}

r3  <- process_tigit("3hMTB")
r24 <- process_tigit("24hMTB")
cat(sprintf("  3h: Low=%d High=%d\n",
            sum(r3$obj$TIGIT_group == "Low"),
            sum(r3$obj$TIGIT_group == "High")))
cat(sprintf("  24h: Low=%d High=%d\n",
            sum(r24$obj$TIGIT_group == "Low"),
            sum(r24$obj$TIGIT_group == "High")))

compute_heat <- function(nk_sub, genes_use) {
  expr_mat <- GetAssayData(nk_sub, layer = "data")
  g_avail <- intersect(genes_use, rownames(expr_mat))
  # Per-gene z-score
  z <- t(scale(t(expr_mat[g_avail, , drop = FALSE])))
  z[z > 3] <- 3
  z[z < -3] <- -3
  # Group means
  lo <- rowMeans(z[, nk_sub$TIGIT_group == "Low", drop = FALSE])
  hi <- rowMeans(z[, nk_sub$TIGIT_group == "High", drop = FALSE])
  m <- cbind(Low = lo, High = hi)
  rownames(m) <- g_avail
  # Wilcoxon test per gene
  wp <- sapply(g_avail, function(g) {
    tryCatch(
      wilcox.test(expr_mat[g, ] ~ nk_sub$TIGIT_group, exact = FALSE)$p.value,
      error = function(e) NA
    )
  })
  list(mat = m, p = wp)
}

h3  <- compute_heat(r3$obj, all_keep)
h24 <- compute_heat(r24$obj, all_keep)

common_genes <- intersect(rownames(h3$mat), rownames(h24$mat))
cat(sprintf("  Common genes: %d\n", length(common_genes)))

# Assign modules
gene_module <- rep("OXPHOS", length(common_genes))
gene_module[common_genes %in% gly_keep] <- "Glycolysis"
gene_module[common_genes %in% lac_keep] <- "Lactate"
names(gene_module) <- common_genes

# Cluster within each module (using 24h data as anchor)
cluster_module <- function(mat, genes_in_module) {
  if (length(genes_in_module) <= 1) return(genes_in_module)
  d <- dist(mat[genes_in_module, , drop = FALSE])
  hc <- hclust(d, method = "ward.D2")
  genes_in_module[hc$order]
}

mod_order <- c("Glycolysis", "Lactate", "OXPHOS")
gene_order <- character()
for (mod in mod_order) {
  g_mod <- names(gene_module)[gene_module == mod]
  if (length(g_mod) > 0) {
    g_clust <- cluster_module(h24$mat, g_mod)
    gene_order <- c(gene_order, g_clust)
  }
}

# Combined heatmap matrix
heat_combined <- cbind(
  h3$mat[gene_order, c("Low", "High")],
  h24$mat[gene_order, c("Low", "High")]
)
colnames(heat_combined) <- c("3h_Low", "3h_High", "24h_Low", "24h_High")
max_abs_c <- max(abs(heat_combined), na.rm = TRUE)

# Module gap positions
mod_present <- unique(gene_module[gene_order])
mod_counts  <- sapply(mod_present, function(m) sum(gene_module[gene_order] == m))
gap_pos <- cumsum(mod_counts[mod_counts > 0])
if (length(gap_pos) > 1) {
  gap_pos <- gap_pos[-length(gap_pos)]
}
# If there is only one module, no gap is needed
if (length(gap_pos) == 0) {
  gap_pos <- integer(0)
}

# Annotations
ann_col_comb <- data.frame(
  Timepoint = c("3h", "3h", "24h", "24h"),
  TIGIT     = c("Low", "High", "Low", "High"),
  row.names = colnames(heat_combined)
)
ann_row_comb <- data.frame(
  Module = gene_module[gene_order],
  row.names = gene_order
)

pdf(file.path(OUT_DIR, "NK_TIGIT_Metabolic_Heatmap.pdf"),
    width = 8, height = 9)
pheatmap(
  heat_combined,
  color             = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  breaks            = seq(-max_abs_c, max_abs_c, length.out = 101),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  annotation_col    = ann_col_comb,
  annotation_row    = ann_row_comb,
  cellwidth         = 28,
  cellheight        = 16,
  gaps_row          = as.numeric(gap_pos),
  gaps_col          = 2,
  annotation_colors = list(
    Timepoint = c("3h" = "#D55E00", "24h" = "#CC79A7"),
    TIGIT     = c(Low = "#E8A838",  High = "#C43A3A"),
    Module    = c(Glycolysis = "#E8A838", Lactate = "#FC4E2A", OXPHOS = "#800026")
  ),
  display_numbers = FALSE,
  fontsize_row    = 10,
  fontsize_col    = 12,
  border_color    = NA,
  main            = "NK: TIGIT Low vs High (3h vs 24h)"
)
dev.off()
cat("  Heatmap saved\n")

# Summary statistics table
stats_df <- data.frame(
  Gene      = gene_order,
  H3_Low    = round(h3$mat[gene_order, "Low"], 3),
  H3_High   = round(h3$mat[gene_order, "High"], 3),
  H24_Low   = round(h24$mat[gene_order, "Low"], 3),
  H24_High  = round(h24$mat[gene_order, "High"], 3),
  H3_Diff   = round(h3$mat[gene_order, "High"] - h3$mat[gene_order, "Low"], 3),
  H24_Diff  = round(h24$mat[gene_order, "High"] - h24$mat[gene_order, "Low"], 3),
  p_adj_24h = round(p.adjust(h24$p[gene_order], "BH"), 4),
  Module    = gene_module[gene_order],
  row.names = NULL
)
stats_df$sig24 <- ifelse(stats_df$p_adj_24h < 0.001, "***",
                  ifelse(stats_df$p_adj_24h < 0.01,  "**",
                  ifelse(stats_df$p_adj_24h < 0.05,  "*",  "ns")))
cat(sprintf("  Significant at 24h (p.adj < 0.05): %d / %d\n",
            sum(stats_df$sig24 != "ns"), nrow(stats_df)))
print(stats_df, row.names = FALSE)
write.csv(stats_df, file.path(OUT_DIR, "NK_TIGIT_Metabolic_Stats.csv"),
          row.names = FALSE)
cat("  Stats CSV saved\n")

# ============================================================
# OUTPUT SUMMARY
# ============================================================
cat("\n=== OUTPUT SUMMARY ===\n")
output_files <- list.files(OUT_DIR, pattern = "\\.pdf$|\\.csv$")
for (f in output_files) {
  cat(sprintf("  %s\n", file.path(OUT_DIR, f)))
}
cat(sprintf("\nTotal: %d file(s) in %s\n", length(output_files), OUT_DIR))
cat("\n=== ALL DONE ===\n")
