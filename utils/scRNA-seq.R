# ==================== utils/scRNA-seq.R ====================
# Shared utility functions for scRNA-seq analysis scripts
# TB-TIGIT-NK manuscript
# ============================================================

#' Compute pathway-level global z-score (reference = all cells at that timepoint)
#'
#' @param nk_full Seurat object containing all cells at a given timepoint
#' @param nk_sub  Seurat object containing TIGIT Low/High subset
#' @param pwy_list Named list of gene vectors, one per pathway
#' @return List with mat (Low/High means of z-scores) and p (Wilcoxon p-values)
compute_pw_global_z <- function(nk_full, nk_sub, pwy_list) {
  expr_full <- GetAssayData(nk_full, layer = "data")
  score_full <- sapply(pwy_list, function(gs) colMeans(expr_full[gs, , drop = FALSE]))
  pop_mean <- colMeans(score_full, na.rm = TRUE)
  pop_sd   <- apply(score_full, 2, sd, na.rm = TRUE)

  expr_sub <- GetAssayData(nk_sub, layer = "data")
  score_sub <- sapply(pwy_list, function(gs) colMeans(expr_sub[gs, , drop = FALSE]))

  z_sub <- sweep(score_sub, 2, pop_mean, "-")
  z_sub <- sweep(z_sub, 2, pop_sd, "/")
  z_sub[z_sub > 3] <- 3; z_sub[z_sub < -3] <- -3

  lo <- colMeans(z_sub[nk_sub$TIGIT_group == "Low", , drop = FALSE])
  hi <- colMeans(z_sub[nk_sub$TIGIT_group == "High", , drop = FALSE])
  m <- cbind(Low = lo, High = hi)

  wp <- sapply(colnames(z_sub), function(pw) {
    tryCatch(wilcox.test(z_sub[, pw] ~ nk_sub$TIGIT_group, exact = FALSE)$p.value,
             error = function(e) NA)
  })
  list(mat = m, p = wp)
}

#' TIGIT grouping: split NK cells into Low/High at a given timepoint
#'
#' @param nk  Seurat object (NK cells)
#' @param tp  Timepoint string (e.g. "3hMTB", "24hMTB")
#' @param n_each Number of cells per group
#' @return Seurat object with TIGIT_group metadata column (Low/High)
process_tigit <- function(nk, tp, n_each = 800) {
  nk_tp <- nk[, nk$timepoint == tp]
  tig <- GetAssayData(nk_tp, layer = "data")["TIGIT", ]
  tig_pos <- tig[tig > 0]
  tig_sorted <- sort(tig_pos)
  if (length(tig_sorted) < n_each * 2) n_each <- floor(length(tig_sorted) / 2)
  low_cells <- names(tig_sorted[1:n_each])
  high_cells <- names(rev(tig_sorted)[1:n_each])
  nk_sub <- nk_tp[, c(low_cells, high_cells)]
  nk_sub$TIGIT_group <- factor(ifelse(colnames(nk_sub) %in% high_cells, "High", "Low"),
                                levels = c("Low", "High"))
  nk_sub
}

#' Build pathway heatmap from TIGIT Low/High comparison (3h + 24h)
#'
#' @param nk          Full NK Seurat object
#' @param pathway_genes Named list of gene vectors
#' @param pathway_cat Named category vector for each pathway
#' @param cat_order   Ordered category names
#' @param out_dir     Output directory for PDF + CSV
#' @param file_prefix Prefix for output filenames
#' @param plot_title  Title for the heatmap
#' @return Invisibly returns the stats data.frame
build_tigit_heatmap <- function(nk, pathway_genes, pathway_cat, cat_order,
                                 out_dir, file_prefix, plot_title) {
  cat("  >> Computing TIGIT Low/High groups...\n")
  r3  <- process_tigit(nk, "3hMTB")
  r24 <- process_tigit(nk, "24hMTB")
  cat(sprintf("  3h: Low=%d High=%d\n", sum(r3$TIGIT_group == "Low"), sum(r3$TIGIT_group == "High")))
  cat(sprintf("  24h: Low=%d High=%d\n", sum(r24$TIGIT_group == "Low"), sum(r24$TIGIT_group == "High")))

  nk_3h_full <- nk[, nk$timepoint == "3hMTB"]
  nk_24h_full <- nk[, nk$timepoint == "24hMTB"]

  cat("  >> Computing global z-scores...\n")
  h3  <- compute_pw_global_z(nk_3h_full, r3, pathway_genes)
  h24 <- compute_pw_global_z(nk_24h_full, r24, pathway_genes)

  common_pw <- intersect(rownames(h3$mat), rownames(h24$mat))
  cat(sprintf("  Common pathways: %d\n", length(common_pw)))

  pw_cat_final <- pathway_cat[common_pw]
  h24_sub <- h24$mat[common_pw, , drop = FALSE]

  gene_order <- character()
  for (ctg in cat_order) {
    pw_c <- names(pw_cat_final)[pw_cat_final == ctg]
    if (length(pw_c) > 1) {
      d <- dist(h24_sub[pw_c, , drop = FALSE])
      pw_c <- pw_c[hclust(d, method = "ward.D2")$order]
    }
    if (length(pw_c) > 0) gene_order <- c(gene_order, pw_c)
  }

  heat_c <- cbind(
    h3$mat[gene_order, c("Low", "High")],
    h24$mat[gene_order, c("Low", "High")]
  )
  colnames(heat_c) <- c("3h_Low", "3h_High", "24h_Low", "24h_High")
  max_abs_c <- max(abs(heat_c), na.rm = TRUE)

  pw_cat_ordered <- pw_cat_final[gene_order]
  cat_counts <- sapply(cat_order, function(x) sum(pw_cat_ordered == x))
  gap_pos <- cumsum(cat_counts[cat_counts > 0])
  if (length(gap_pos) > 1) gap_pos <- gap_pos[-length(gap_pos)]

  ann_col <- data.frame(
    Timepoint = c("3h", "3h", "24h", "24h"),
    TIGIT = c("Low", "High", "Low", "High"),
    row.names = colnames(heat_c)
  )
  ann_row <- data.frame(Category = factor(pw_cat_ordered, levels = cat_order),
                        row.names = gene_order)

  cat_colors <- c(
    "Energy Metabolism" = "#E64B35", "Amino Acid" = "#4DBBD5",
    "Lipid Metabolism" = "#F39B7F", "Nucleotide" = "#00A087",
    "Redox" = "#3C5488"
  )

  cat("  >> Plotting heatmap...\n")
  pdf(file.path(out_dir, paste0(file_prefix, ".pdf")), width = 7, height = max(4, length(gene_order) * 0.35 + 2))
  pheatmap(heat_c,
    color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    breaks = seq(-max_abs_c, max_abs_c, length.out = 101),
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_col = ann_col, annotation_row = ann_row,
    cellwidth = 24, cellheight = max(14, min(20, 200 / length(gene_order))),
    gaps_row = as.numeric(gap_pos), gaps_col = 2,
    annotation_colors = list(
      Timepoint = c("3h" = "#D55E00", "24h" = "#CC79A7"),
      TIGIT = c(Low = "#E8A838", High = "#C43A3A"),
      Category = cat_colors
    ),
    display_numbers = FALSE, fontsize_row = 10, fontsize_col = 12, border_color = NA,
    main = plot_title
  )
  dev.off()
  cat("  >> Saved: ", file.path(out_dir, paste0(file_prefix, ".pdf")), "\n")

  stats_df <- data.frame(
    Pathway = gene_order,
    H3_Low = round(heat_c[gene_order, "3h_Low"], 3),
    H3_High = round(heat_c[gene_order, "3h_High"], 3),
    H24_Low = round(heat_c[gene_order, "24h_Low"], 3),
    H24_High = round(heat_c[gene_order, "24h_High"], 3),
    Diff_3h = round(heat_c[gene_order, "3h_High"] - heat_c[gene_order, "3h_Low"], 3),
    Diff_24h = round(heat_c[gene_order, "24h_High"] - heat_c[gene_order, "24h_Low"], 3),
    Category = pw_cat_ordered, row.names = NULL
  )
  stats_df$p_adj_3h  <- round(p.adjust(h3$p[gene_order], "BH"), 4)
  stats_df$p_adj_24h <- round(p.adjust(h24$p[gene_order], "BH"), 4)
  stats_df$sig <- ifelse(stats_df$p_adj_24h < 0.001, "***",
                  ifelse(stats_df$p_adj_24h < 0.01, "**",
                  ifelse(stats_df$p_adj_24h < 0.05, "*", "ns")))

  write.csv(stats_df, file.path(out_dir, paste0(file_prefix, "_Stats.csv")), row.names = FALSE)
  cat("  >> Saved stats: ", file.path(out_dir, paste0(file_prefix, "_Stats.csv")), "\n")

  print(stats_df, row.names = FALSE)
  invisible(stats_df)
}
