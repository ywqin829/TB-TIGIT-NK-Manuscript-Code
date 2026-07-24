# =============================================================================
# GSVA_1.R -- GSVA analysis (KEGG + GO) of bulk RNA-seq data
# TB-TIGIT-NK manuscript
#
# Input:  ./data/exp.csv  (normalized expression matrix, genes x samples)
# Output: ./output/
#   1. gsva_kegg_result.csv             — GSVA KEGG scores
#   2. gsva_go_result.csv               — GSVA GO scores
#   3. GSVA_KEGG_volcano.pdf            — Volcano plot (KEGG)
#   4. GSVA_KEGG_heatmap.pdf            — Heatmap of top KEGG pathways
#   5. GSVA_KEGG_barplot_v1.pdf         — Barplot: significant KEGG pathways only
#   6. GSVA_KEGG_barplot.pdf            — Barplot: top KEGG by logFC (sig colored)
#   7. GSVA_GO_volcano.pdf              — Volcano plot (GO)
#   8. GSVA_GO_heatmap.pdf              — Heatmap of top GO pathways
#   9. GSVA_GO_barplot.pdf              — Barplot: top GO by logFC (sig colored)
#
# Notes:
#   - GSVA computation is wrapped in if (FALSE) to skip after first run
#   - On subsequent runs, pre-computed CSV files are loaded
#   - Sample group assignment: "Ctrl*" -> IgG, "TIGIT*" -> alphaTIGIT
# =============================================================================

# ===== SETUP =====
rm(list = ls())
cat("=== GSVA Analysis ===\n")

# -- Libraries --
library(msigdbr)
library(dplyr)
library(GSVA)
library(BiocParallel)
library(limma)
library(ggplot2)
library(pheatmap)
library(tibble)
library(tidyr)
library(ggthemes)
library(ggprism)
library(grid)

# -- Directories --
OUT_DIR <- "./output"
DATA_DIR <- "./data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===== LOAD DATA =====
cat("\n===== Loading expression data =====\n")
exp <- read.csv(file.path(DATA_DIR, "exp.csv"), row.names = 1)
cat(sprintf("  Genes: %d, Samples: %d\n", nrow(exp), ncol(exp)))

# ===== Gene set preparation function =====
prepare_gene_sets <- function(species, collection, subcollection = NULL) {
  if (!is.null(subcollection)) {
    df_all <- do.call(rbind, lapply(subcollection, function(subcol) {
      msigdbr(species = species, collection = collection, subcollection = subcol)
    }))
  } else {
    df_all <- msigdbr(species = species, collection = collection)
  }
  gene_df <- df_all %>% dplyr::select(gs_name, gene_symbol)
  gene_list <- split(gene_df$gene_symbol, gene_df$gs_name)
  gene_list <- gene_list[sapply(gene_list, length) >= 10]
  return(gene_list)
}

# ===== Prepare group info (shared across KEGG and GO) =====
samples_all <- colnames(exp)
group_info <- data.frame(
  Sample = samples_all,
  Group = ifelse(grepl("^Ctrl", samples_all), "IgG",
                 ifelse(grepl("^TIGIT", samples_all), "alphaTIGIT", NA))
)
group_info <- subset(group_info, Group %in% c("IgG", "alphaTIGIT"))
group_info <- group_info[order(factor(group_info$Group, levels = c("IgG", "alphaTIGIT"))), ]
cat(sprintf("  Samples: %d (IgG: %d, alphaTIGIT: %d)\n",
            nrow(group_info),
            sum(group_info$Group == "IgG"),
            sum(group_info$Group == "alphaTIGIT")))

# =============================================================================
# PART 1: GSVA KEGG Computation (expensive, one-time)
# =============================================================================
cat("\n===== PART 1: GSVA KEGG =====\n")

if (FALSE) {
  cat("  Preparing KEGG gene sets...\n")
  kegg_list <- prepare_gene_sets(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:KEGG_LEGACY"
  )
  cat(sprintf("  KEGG pathways (>=10 genes): %d\n", length(kegg_list)))

  dat <- as.matrix(exp)
  rownames(dat) <- toupper(rownames(dat))
  kegg_list <- lapply(kegg_list, toupper)

  register(MulticoreParam(workers = 4))
  param <- gsvaParam(exprData = dat, geneSets = kegg_list, kcdf = "Gaussian")
  cat("  Running GSVA KEGG (this may take a while)...\n")
  gsva_kegg_result <- gsva(param, verbose = TRUE,
                           BPPARAM = MulticoreParam(workers = 4))

  write.csv(gsva_kegg_result, file.path(OUT_DIR, "gsva_kegg_result.csv"))
  cat("  Saved: ", file.path(OUT_DIR, "gsva_kegg_result.csv"), "\n")
} else {
  cat("  Loading pre-computed GSVA KEGG scores...\n")
  gsva_kegg_result <- read.csv(file.path(OUT_DIR, "gsva_kegg_result.csv"),
                                row.names = 1)
  cat(sprintf("  Pathways: %d, Samples: %d\n", nrow(gsva_kegg_result),
              ncol(gsva_kegg_result)))
}

# =============================================================================
# PART 2: GSVA KEGG Differential Analysis
# =============================================================================
cat("\n===== PART 2: KEGG Diff Analysis =====\n")

# -- Reorder columns to match group_info --
gsva_kegg_result <- gsva_kegg_result[, group_info$Sample]

group_list <- factor(group_info$Group, levels = c("IgG", "alphaTIGIT"))

# -- Design matrix --
design <- model.matrix(~ group_list)

# -- limma differential analysis --
fit <- lmFit(gsva_kegg_result, design)
fit <- eBayes(fit)
diff_gsva_kegg_result <- as.data.frame(
  topTable(fit, coef = 2, number = Inf, adjust = "fdr")
)

# -- Thresholds --
LogFC <- 0.5
adj.P.Val <- 0.05

diff_gsva_kegg_result$change <- ifelse(
  (diff_gsva_kegg_result$adj.P.Val < adj.P.Val) &
    (diff_gsva_kegg_result$logFC < -LogFC), "Down",
  ifelse((diff_gsva_kegg_result$adj.P.Val < adj.P.Val) &
           (diff_gsva_kegg_result$logFC > LogFC), "Up", "NS")
)

cat("  Change distribution:\n")
print(table(diff_gsva_kegg_result$change))

# ---- Volcano plot (KEGG) ----
cat("  Plotting KEGG volcano...\n")
p <- ggplot(data = diff_gsva_kegg_result,
            aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.6, size = 2.0, aes(color = change)) +
  ylab("-Log10(P.value)") +
  scale_color_manual(
    values = c("Down" = "#0000EE", "NS" = "grey", "Up" = "#8B2323")
  ) +
  geom_vline(xintercept = c(-LogFC, LogFC), lty = 4,
             col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4,
             col = "black", lwd = 0.8) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 14))
ggsave(file.path(OUT_DIR, "GSVA_KEGG_volcano.pdf"), plot = p,
       width = 8, height = 6)
cat("  Saved: ", file.path(OUT_DIR, "GSVA_KEGG_volcano.pdf"), "\n")

# ---- Heatmap (KEGG) ----
cat("  Plotting KEGG heatmap...\n")
sig_gsva_kegg_result <- diff_gsva_kegg_result[
  diff_gsva_kegg_result$adj.P.Val < 0.05 &
    abs(diff_gsva_kegg_result$logFC) > 0.5, ]
sig_gsva_kegg_result <- sig_gsva_kegg_result[
  order(sig_gsva_kegg_result$adj.P.Val,
        -abs(sig_gsva_kegg_result$logFC)), ]

top_n <- 20
top_sig_gsva_kegg_result <- head(rownames(sig_gsva_kegg_result), top_n)

heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
heatmap_matrix <- gsva_kegg_result[top_sig_gsva_kegg_result, ]

annotation_col <- data.frame(Group = group_list)
rownames(annotation_col) <- group_info$Sample

pdf(file.path(OUT_DIR, "GSVA_KEGG_heatmap.pdf"), width = 10, height = 8)
pheatmap(
  heatmap_matrix,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cellwidth = 30,
  cellheight = 25,
  scale = "row",
  main = "Heatmap of Top GSVA Results (KEGG)",
  labels_row = gsub("KEGG_", "", top_sig_gsva_kegg_result),
  annotation_col = annotation_col,
  fontsize = 19,
  fontsize_row = 17
)
dev.off()
cat("  Saved: ", file.path(OUT_DIR, "GSVA_KEGG_heatmap.pdf"), "\n")

# ---- Barplot (KEGG), Version 1: Significant pathways only ----
cat("  Plotting KEGG barplot (significant only)...\n")
deg <- diff_gsva_kegg_result

upregulated <- deg %>%
  filter(logFC > 0.5, adj.P.Val < 0.05) %>%
  arrange(desc(logFC)) %>%
  dplyr::slice(1:min(25, n()))
downregulated <- deg %>%
  filter(logFC < 0, adj.P.Val < 0.05) %>%
  arrange(logFC) %>%
  dplyr::slice(1:min(23, n()))
Diff <- rbind(upregulated, downregulated)

dat_plot <- data.frame(
  id = row.names(Diff),
  p = Diff$adj.P.Val,
  lgfc = Diff$logFC
)
dat_plot$group <- ifelse(dat_plot$lgfc > 0, 1, -1)
dat_plot$lg_p <- -log10(dat_plot$p) * dat_plot$group
dat_plot$id <- gsub("KEGG_", "", dat_plot$id)
dat_plot$threshold <- factor(
  ifelse(dat_plot$p <= 0.05,
         ifelse(dat_plot$lgfc > 0, "Up", "Down"), "Not"),
  levels = c("Up", "Down", "Not")
)
dat_plot <- dat_plot[order(dat_plot$lg_p), ]
dat_plot$id <- factor(dat_plot$id, levels = dat_plot$id)

p <- ggplot(data = dat_plot, aes(x = id, y = lg_p, fill = threshold)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("Up" = "#FF0000", "Down" = "#0000FF", "Not" = "#cccccc")
  ) +
  geom_hline(yintercept = c(-1.3, 1.3), color = "white",
             linewidth = 0.5, lty = "dashed") +
  xlab("") +
  ylab("-log10(adj.P.Value) of GSVA KEGG score") +
  guides(fill = "none") +
  theme_prism(border = TRUE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14)
  ) +
  geom_text(
    aes(
      y = ifelse(lg_p > 0, -0.02 * max(abs(lg_p)),
                 max(abs(lg_p)) * 0.02),
      label = id,
      hjust = ifelse(lg_p > 0, 1, 0)
    ),
    color = ifelse(dat_plot$p <= 0.05, "black", "grey"),
    size = 3, vjust = 0.5
  )
ggsave(file.path(OUT_DIR, "GSVA_KEGG_barplot_v1.pdf"), plot = p,
       width = 10, height = 8)
cat("  Saved: ", file.path(OUT_DIR, "GSVA_KEGG_barplot_v1.pdf"), "\n")

# ---- Barplot (KEGG), Version 2: Top by logFC (significance colored) ----
cat("  Plotting KEGG barplot (top by logFC)...\n")
p_cutoff <- 0.05

upregulated <- deg %>%
  filter(logFC > 0) %>%
  arrange(desc(logFC)) %>%
  slice(1:min(25, n()))
downregulated <- deg %>%
  filter(logFC < 0) %>%
  arrange(logFC) %>%
  slice(1:min(25, n()))
Diff <- rbind(upregulated, downregulated)

dat_plot <- data.frame(
  id = row.names(Diff),
  p = Diff$adj.P.Val,
  lgfc = Diff$logFC
)
dat_plot$group <- ifelse(dat_plot$lgfc > 0, 1, -1)
dat_plot$lg_p <- -log10(dat_plot$p) * dat_plot$group
dat_plot$id <- gsub("KEGG_", "", dat_plot$id)
dat_plot$threshold <- factor(
  ifelse(dat_plot$p < p_cutoff,
         ifelse(dat_plot$lgfc > 0, "Up", "Down"), "Not"),
  levels = c("Up", "Down", "Not")
)
dat_plot <- dat_plot[order(dat_plot$lg_p), ]
dat_plot$id <- factor(dat_plot$id, levels = dat_plot$id)

p <- ggplot(data = dat_plot, aes(x = id, y = lg_p, fill = threshold)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c("Up" = "#FF0000", "Down" = "#0000FF", "Not" = "#CCCCCC")
  ) +
  geom_hline(yintercept = c(-1.3, 1.3), color = "white",
             linewidth = 0.5, lty = "dashed") +
  xlab("") +
  ylab("-log10(adj.P.Value) of GSVA KEGG score") +
  guides(fill = "none") +
  theme_prism(border = TRUE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)
  ) +
  geom_text(
    aes(
      y = ifelse(lg_p > 0, -0.02 * max(abs(lg_p)),
                 max(abs(lg_p)) * 0.02),
      label = id,
      hjust = ifelse(lg_p > 0, 1, 0)
    ),
    color = ifelse(dat_plot$p < p_cutoff, "black", "grey"),
    size = 4, vjust = 0.5
  )
ggsave(file.path(OUT_DIR, "GSVA_KEGG_barplot.pdf"), plot = p,
       width = 10, height = 8)
cat("  Saved: ", file.path(OUT_DIR, "GSVA_KEGG_barplot.pdf"), "\n")

# =============================================================================
# PART 3: GSVA GO Computation (expensive, one-time)
# =============================================================================
cat("\n===== PART 3: GSVA GO =====\n")

if (FALSE) {
  cat("  Preparing GO gene sets...\n")
  go_list <- prepare_gene_sets(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = c("GO:BP", "GO:CC", "GO:MF")
  )
  cat(sprintf("  GO terms (>=10 genes): %d\n", length(go_list)))

  dat <- as.matrix(exp)
  rownames(dat) <- toupper(rownames(dat))
  go_list <- lapply(go_list, toupper)

  register(MulticoreParam(workers = 4))
  param_go <- gsvaParam(exprData = dat, geneSets = go_list, kcdf = "Gaussian")
  cat("  Running GSVA GO (this may take a while)...\n")
  gsva_go_result <- gsva(param_go, verbose = TRUE,
                         BPPARAM = MulticoreParam(workers = 4))

  write.csv(gsva_go_result, file.path(OUT_DIR, "gsva_go_result.csv"))
  cat("  Saved: ", file.path(OUT_DIR, "gsva_go_result.csv"), "\n")
} else {
  cat("  Loading pre-computed GSVA GO scores...\n")
  gsva_go_result <- read.csv(file.path(OUT_DIR, "gsva_go_result.csv"),
                              row.names = 1)
  cat(sprintf("  Pathways: %d, Samples: %d\n", nrow(gsva_go_result),
              ncol(gsva_go_result)))
}

# =============================================================================
# PART 4: GSVA GO Differential Analysis
# =============================================================================
cat("\n===== PART 4: GO Diff Analysis =====\n")

# -- Reorder columns --
gsva_go_result <- gsva_go_result[, group_info$Sample]

group_list <- factor(group_info$Group, levels = c("IgG", "alphaTIGIT"))

# -- Design matrix --
design <- model.matrix(~ group_list)

# -- limma differential analysis --
fit <- lmFit(gsva_go_result, design)
fit <- eBayes(fit)
diff_gsva_go_result <- as.data.frame(
  topTable(fit, coef = 2, number = Inf, adjust = "fdr")
)

# -- Thresholds (GO uses |logFC| > 1, stricter than KEGG) --
logFC <- 1
adj.P.Val <- 0.05

diff_gsva_go_result$change <- ifelse(
  (diff_gsva_go_result$adj.P.Val < adj.P.Val) &
    (diff_gsva_go_result$logFC < -logFC), "Down",
  ifelse((diff_gsva_go_result$adj.P.Val < adj.P.Val) &
           (diff_gsva_go_result$logFC > logFC), "Up", "NS")
)

cat("  Change distribution:\n")
print(table(diff_gsva_go_result$change))

# ---- Volcano plot (GO) ----
cat("  Plotting GO volcano...\n")
p <- ggplot(data = diff_gsva_go_result,
            aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.6, size = 2.0, aes(color = change)) +
  ylab("-Log10(adj.P.Val)") +
  scale_color_manual(
    values = c("Down" = "#0000EE", "NS" = "grey", "Up" = "#8B2323")
  ) +
  geom_vline(xintercept = c(-logFC, logFC), lty = 4,
             col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(adj.P.Val), lty = 4,
             col = "black", lwd = 0.8) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 14))
ggsave(file.path(OUT_DIR, "GSVA_GO_volcano.pdf"), plot = p,
       width = 8, height = 6)
cat("  Saved: ", file.path(OUT_DIR, "GSVA_GO_volcano.pdf"), "\n")

# ---- Heatmap (GO) ----
cat("  Plotting GO heatmap...\n")
sig_go_pathways <- diff_gsva_go_result %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)
top_go <- sig_go_pathways %>%
  arrange(logFC) %>%
  dplyr::slice(1:20)

heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
heatmap_matrix <- gsva_go_result[rownames(top_go), ]

annotation_col <- data.frame(Group = group_list)
rownames(annotation_col) <- group_info$Sample

pdf(file.path(OUT_DIR, "GSVA_GO_heatmap.pdf"), width = 10, height = 8)
pheatmap(
  mat = heatmap_matrix,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cellwidth = 20,
  cellheight = 20,
  fontsize = 10,
  fontsize_row = 14,
  scale = "row",
  main = "Heatmap of GSVA GO (adj.P.Val < 0.05, |logFC| > 1.0)",
  annotation_col = annotation_col,
  labels_row = gsub("^GO[_:]*", "", rownames(heatmap_matrix))
)
dev.off()
cat("  Saved: ", file.path(OUT_DIR, "GSVA_GO_heatmap.pdf"), "\n")

# ---- Barplot (GO): Top by logFC (significance colored) ----
# Note: uses diff_gsva_go_result (GO results, not KEGG)
cat("  Plotting GO barplot...\n")
p_cutoff <- 0.05
logFC_cutoff <- 1.0
n_top <- 20

selected_paths <- diff_gsva_go_result %>%
  dplyr::filter(abs(logFC) > logFC_cutoff & adj.P.Val < p_cutoff) %>%
  mutate(id = row.names(.)) %>%
  group_by(sign(logFC)) %>%
  arrange(desc(abs(logFC))) %>%
  dplyr::slice(1:min(n_top, n())) %>%
  ungroup()

cat("  Selected GO paths (head):\n")
print(head(selected_paths, 10))
cat(sprintf("  Up: %d, Down: %d\n",
            sum(selected_paths$logFC > 0),
            sum(selected_paths$logFC < 0)))

dat_plot <- data.frame(
  id = selected_paths$id,
  p = selected_paths$adj.P.Val,
  lgfc = selected_paths$logFC
)
dat_plot$group <- ifelse(dat_plot$lgfc > 0, 1, -1)
dat_plot$lg_p <- -log10(dat_plot$p) * dat_plot$group
dat_plot$id <- gsub("GO_", "", dat_plot$id)
dat_plot$threshold <- factor(
  ifelse(dat_plot$p < p_cutoff,
         ifelse(dat_plot$lgfc > 0, "Up", "Down"), "Not"),
  levels = c("Up", "Down", "Not")
)

dat_plot <- dat_plot[order(dat_plot$lg_p), ]
dat_plot$id <- factor(dat_plot$id, levels = dat_plot$id)

p <- ggplot(data = dat_plot, aes(x = id, y = lg_p, fill = threshold)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c("Up" = "#FF0000", "Down" = "#0000FF", "Not" = "#CCCCCC")
  ) +
  geom_hline(yintercept = c(-log10(p_cutoff), log10(p_cutoff)),
             color = "white", linewidth = 0.5, lty = "dashed") +
  xlab("") +
  ylab("-log10(adj.P.Value) of GSVA_GO score") +
  guides(fill = "none") +
  theme_prism(border = TRUE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)
  ) +
  geom_text(
    aes(
      y = ifelse(lg_p > 0, -0.03 * max(abs(lg_p)),
                 max(abs(lg_p)) * 0.03),
      label = id,
      hjust = ifelse(lg_p > 0, 1, 0)
    ),
    color = ifelse(dat_plot$p < p_cutoff, "black", "grey"),
    size = 3, vjust = 0.5
  )
ggsave(file.path(OUT_DIR, "GSVA_GO_barplot.pdf"), plot = p,
       width = 10, height = 8)
cat("  Saved: ", file.path(OUT_DIR, "GSVA_GO_barplot.pdf"), "\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  GSVA Analysis Complete\n")
cat("================================================================\n")
cat("  Output files:\n")
cat("    ", file.path(OUT_DIR, "gsva_kegg_result.csv"), "\n")
cat("    ", file.path(OUT_DIR, "gsva_go_result.csv"), "\n")
cat("    ", file.path(OUT_DIR, "GSVA_KEGG_volcano.pdf"), "\n")
cat("    ", file.path(OUT_DIR, "GSVA_KEGG_heatmap.pdf"), "\n")
cat("    ", file.path(OUT_DIR, "GSVA_KEGG_barplot_v1.pdf"), "\n")
cat("    ", file.path(OUT_DIR, "GSVA_KEGG_barplot.pdf"), "\n")
cat("    ", file.path(OUT_DIR, "GSVA_GO_volcano.pdf"), "\n")
cat("    ", file.path(OUT_DIR, "GSVA_GO_heatmap.pdf"), "\n")
cat("    ", file.path(OUT_DIR, "GSVA_GO_barplot.pdf"), "\n")
cat("================================================================\n")
