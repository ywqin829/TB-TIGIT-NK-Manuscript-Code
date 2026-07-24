# ==============================================================================
# TIGIT代谢组学测序分析.R — Metabolomics analysis: TIGIT High vs Low
# TB-TIGIT-NK manuscript
#
# Input:  111.csv (metabolomics data matrix, rows = metabolites, cols = samples)
# Output panels:
#   1. Normalization + DE (limma + VIP scores via ropls)
#   2. Volcano plot of differential metabolites
#   3. Heatmap of top differential metabolites
#   4. KEGG pathway enrichment (MSEA)
#   5. GSEA analysis
#   6. CIBERSORT immune infiltration
# =============================================================================

# ===== SETUP =====
rm(list = ls())
cat("=== TIGIT Metabolomics Analysis ===\n")

library(limma)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(RColorBrewer)
library(FactoMineR)
library(factoextra)
library(ropls)
library(KEGGREST)
library(reshape2)
library(ggsignif)

source("../utils/go_kegg.R")

OUT_DIR <- "./output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
DATA_DIR <- "./data"

set.seed(42)

# ===== LOAD DATA =====
cat("\n===== Loading metabolomics data =====\n")
metab <- read.csv(file.path(DATA_DIR, "111.csv"), header = TRUE, row.names = 1)
cat(sprintf("  Raw dimensions: %d metabolites x %d samples\n", nrow(metab), ncol(metab)))

# Filter low-abundance metabolites (at least 2 non-zero samples)
metab <- metab[rowSums(metab != 0) > 1, ]
cat(sprintf("  After filtering low-count: %d metabolites\n", nrow(metab)))

# ===== SAMPLE GROUPING =====
cat("\n===== Sample grouping =====\n")
samples <- colnames(metab)
group_info <- data.frame(
  Sample = samples,
  Group = ifelse(grepl("^IgG|^Ctrl", samples), "IgG",
         ifelse(grepl("^TIGIT", samples), "alphaTIGIT", NA)),
  stringsAsFactors = FALSE
)
group_info <- subset(group_info, !is.na(Group))
group_info <- group_info[order(group_info$Group), ]
group_list <- factor(group_info$Group, levels = c("IgG", "alphaTIGIT"))

# Subset metabolomics matrix to matched samples
metab <- metab[, group_info$Sample]
cat(sprintf("  IgG: %d, alphaTIGIT: %d\n",
    sum(group_list == "IgG"), sum(group_list == "alphaTIGIT")))

# ===== PART 1: Normalization and Differential Expression =====
# (Expensive — wrapped in if(FALSE), loads saved CSV on re-run)
cat("\n===== PART 1: Differential Expression + VIP Scores =====\n")

if (FALSE) {
  # --- 1a. Quantile normalization ---
  metab_norm <- normalizeBetweenArrays(as.matrix(metab), method = "quantile")
  cat("  Quantile normalization applied.\n")

  # --- 1b. Fold change + t-test ---
  res_df <- apply(metab_norm, 1, function(x) {
    x1 <- x[group_info$Group == "alphaTIGIT"]
    x2 <- x[group_info$Group == "IgG"]
    fc  <- mean(x1) / mean(x2)
    l2  <- log2(fc)
    pv  <- suppressWarnings(t.test(x1, x2)$p.value)
    c(FC = fc, log2FC = l2, P.Value = pv)
  })
  res_df <- as.data.frame(t(res_df))
  res_df$Metabolite <- rownames(res_df)
  res_df$FDR <- p.adjust(res_df$P.Value, method = "fdr")

  # --- 1c. PLS-DA VIP scores via ropls ---
  pls_model <- opls(t(metab_norm), group_info$Group,
                    predI = 1, orthoI = 1, crossvalI = 6)
  vip_scores <- getVipVn(pls_model)
  res_df$VIP <- vip_scores[rownames(res_df)]

  # --- 1d. Classification by VIP + P.Value (original criteria) ---
  res_df$type <- "insig"
  res_df$type[res_df$VIP > 1 & res_df$P.Value < 0.05 & res_df$log2FC > 0] <- "up"
  res_df$type[res_df$VIP > 1 & res_df$P.Value < 0.05 & res_df$log2FC < 0] <- "down"

  # Reorder columns
  res_df <- res_df[, c("Metabolite", "VIP", "P.Value", "FDR", "FC", "log2FC", "type")]

  write.csv(res_df, file.path(OUT_DIR, "metabolomics_DE_results.csv"), row.names = FALSE)
  cat("  Saved: ", file.path(OUT_DIR, "metabolomics_DE_results.csv"), "\n")
} else {
  res_df <- read.csv(file.path(OUT_DIR, "metabolomics_DE_results.csv"),
                     header = TRUE, stringsAsFactors = FALSE)
  # Re-derive type column from stored values
  res_df$type <- "insig"
  res_df$type[res_df$VIP > 1 & res_df$P.Value < 0.05 & res_df$log2FC > 0] <- "up"
  res_df$type[res_df$VIP > 1 & res_df$P.Value < 0.05 & res_df$log2FC < 0] <- "down"
}
cat(sprintf("  Up: %d, Down: %d, Insig: %d\n",
    sum(res_df$type == "up"), sum(res_df$type == "down"), sum(res_df$type == "insig")))

# ===== INTERLUDE: PCA Plot =====
cat("\n===== PCA of metabolomics =====\n")
pca_res <- prcomp(t(metab), scale. = TRUE)
explained_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
p <- fviz_pca_ind(pca_res,
             geom.ind = "point",
             col.ind = group_info$Group,
             palette = c("#4DBBD5", "#E64B35"),
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Group") +
  labs(title = "PCA of Metabolomics (IgG vs alphaTIGIT)",
       x = paste0("PC1 (", explained_var[1], "%)"),
       y = paste0("PC2 (", explained_var[2], "%)"))
ggsave(file.path(OUT_DIR, "metabolomics_PCA.pdf"), p, width = 7, height = 6)
cat("  Saved: ", file.path(OUT_DIR, "metabolomics_PCA.pdf"), "\n")

# ===== PART 2: Volcano Plot =====
cat("\n===== PART 2: Volcano Plot =====\n")

# Prepare key-value for labelled metabolites (top 10 up + top 10 down)
sig_up   <- res_df %>% filter(type == "up")   %>% arrange(P.Value) %>% slice_head(n = 10)
sig_down <- res_df %>% filter(type == "down") %>% arrange(P.Value) %>% slice_head(n = 10)
label_df <- bind_rows(sig_up, sig_down)

p <- EnhancedVolcano(res_df,
                lab            = res_df$Metabolite,
                x              = "log2FC",
                y              = "P.Value",
                selectLab      = label_df$Metabolite,
                pCutoff        = 0.05,
                FCcutoff       = 1,
                pointSize      = c(ifelse(res_df$type != "insig", 3, 1)),
                labSize        = 4.5,
                colAlpha       = 0.8,
                colCustom      = c("insig" = "grey70", "up" = "#E64B35", "down" = "#4DBBD5"),
                legendPosition = "top",
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                title          = "Metabolomics: TIGIT High vs Low",
                subtitle       = 'VIP > 1 & P.Value < 0.05',
                caption        = paste0("Total metabolites: ", nrow(res_df)))
ggsave(file.path(OUT_DIR, "metabolomics_volcano.pdf"), p, width = 10, height = 8)
cat("  Saved: ", file.path(OUT_DIR, "metabolomics_volcano.pdf"), "\n")

# ===== PART 3: Heatmap =====
cat("\n===== PART 3: Heatmap of significant metabolites =====\n")

sig_metabolites <- res_df$Metabolite[res_df$type != "insig"]
if (length(sig_metabolites) > 0) {
  # Expression matrix for significant metabolites
  exp_sig <- metab[sig_metabolites, group_info$Sample, drop = FALSE]

  # Annotation column
  annotation_col <- data.frame(Group = group_list, row.names = group_info$Sample)
  ann_colors <- list(Group = c(IgG = "#4DBBD5", alphaTIGIT = "#E64B35"))

  # --- 3a. Heatmap of all significant metabolites ---
  if (length(sig_metabolites) > 1) {
    pheatmap(as.matrix(exp_sig),
             scale = "row",
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             show_colnames = FALSE,
             clustering_method = "ward.D2",
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
             cellwidth = 25,
             cellheight = 10,
             fontsize = 12,
             main = "Significant Metabolites (VIP > 1, P < 0.05)",
             filename = file.path(OUT_DIR, "metabolomics_heatmap_sig.pdf"))
    cat("  Saved: ", file.path(OUT_DIR, "metabolomics_heatmap_sig.pdf"), "\n")
  }

  # --- 3b. Heatmap: top 50 by absolute log2FC ---
  top50 <- res_df %>%
    arrange(desc(abs(log2FC))) %>%
    slice_head(n = 50) %>%
    pull(Metabolite)
  exp_top <- metab[top50, group_info$Sample, drop = FALSE]
  pheatmap(as.matrix(exp_top),
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_colnames = FALSE,
           clustering_method = "ward.D2",
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           cellwidth = 15,
           cellheight = 10,
           fontsize_row = 9,
           main = "Top 50 Metabolites by |log2FC|",
           filename = file.path(OUT_DIR, "metabolomics_heatmap_top50.pdf"))
  cat("  Saved: ", file.path(OUT_DIR, "metabolomics_heatmap_top50.pdf"), "\n")

  # --- 3c. Boxplots of top significant metabolites ---
  cat("\n  Generating boxplots of significant metabolites...\n")
  exp_long <- exp_sig %>%
    as.data.frame() %>%
    rownames_to_column("Metabolite") %>%
    pivot_longer(-Metabolite, names_to = "Sample", values_to = "Concentration") %>%
    left_join(group_info, by = "Sample")

  p <- ggplot(exp_long, aes(x = Group, y = Concentration, fill = Group)) +
    stat_summary(fun = mean, geom = "bar", width = 0.6, alpha = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
    facet_wrap(~ Metabolite, scales = "free_y", ncol = 5) +
    labs(x = NULL, y = "Concentration",
         title = "Significant Metabolites: IgG vs alphaTIGIT") +
    scale_fill_manual(values = c(IgG = "#4DBBD5", alphaTIGIT = "#E64B35")) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9),
          axis.title = element_text(size = 11),
          legend.position = "none") +
    geom_signif(comparisons = list(c("IgG", "alphaTIGIT")),
                test = "t.test", map_signif_level = FALSE,
                textsize = 3, tip_length = 0.02)
  ggsave(file.path(OUT_DIR, "metabolomics_boxplots.pdf"), p,
         width = 16, height = ceiling(length(sig_metabolites) / 5) * 3.5)
  cat("  Saved: ", file.path(OUT_DIR, "metabolomics_boxplots.pdf"), "\n")
} else {
  cat("  No significant metabolites to plot.\n")
}

# ===== PART 4: KEGG Pathway Enrichment (MSEA) =====
cat("\n===== PART 4: KEGG Pathway Enrichment (MSEA) =====\n")

# Map metabolite names to KEGG compound IDs via KEGGREST
get_kegg_id <- function(metabolite_name) {
  tryCatch({
    res <- keggFind("compound", metabolite_name)
    if (length(res) > 0) return(sub("^cpd:", "", names(res)[1]))
    return(NA)
  }, error = function(e) return(NA))
}

if (FALSE) {
  sig_df <- res_df %>% filter(type != "insig")
  sig_df$KEGG_ID <- sapply(sig_df$Metabolite, get_kegg_id)
  write.csv(sig_df, file.path(OUT_DIR, "metabolomics_sig_with_kegg.csv"), row.names = FALSE)
  cat("  Saved: ", file.path(OUT_DIR, "metabolomics_sig_with_kegg.csv"), "\n")
} else {
  kegg_map_path <- file.path(OUT_DIR, "metabolomics_sig_with_kegg.csv")
  if (file.exists(kegg_map_path)) {
    sig_df <- read.csv(kegg_map_path, stringsAsFactors = FALSE)
  } else {
    sig_df <- res_df %>% filter(type != "insig")
    sig_df$KEGG_ID <- sapply(sig_df$Metabolite, get_kegg_id)
    write.csv(sig_df, kegg_map_path, row.names = FALSE)
  }
}

# Manual MSEA: fetch KEGG pathways and do hypergeometric test
DE_kegg <- sig_df$KEGG_ID[!is.na(sig_df$KEGG_ID)] %>% unique()
cat(sprintf("  Metabolites with valid KEGG IDs: %d\n", length(DE_kegg)))

if (length(DE_kegg) >= 3) {
  # Fetch all human KEGG pathway IDs
  pathway_list <- keggList("pathway", "hsa")
  pathway_ids  <- names(pathway_list)

  # Get compounds for each pathway
  get_compounds <- function(pid) {
    tryCatch({
      detail <- keggGet(pid)[[1]]
      if (!is.null(detail$COMPOUND)) return(names(detail$COMPOUND))
      return(NULL)
    }, error = function(e) NULL)
  }
  pw_compound_list <- lapply(pathway_ids, get_compounds)
  names(pw_compound_list) <- pathway_ids
  pw_compound_list <- pw_compound_list[!sapply(pw_compound_list, is.null)]

  # Background universe
  all_compounds <- unique(unlist(pw_compound_list))
  bg_n <- length(all_compounds)
  de_n <- length(DE_kegg)

  # Hypergeometric test
  msea_list <- lapply(names(pw_compound_list), function(pid) {
    pw_cmpds <- pw_compound_list[[pid]]
    overlap  <- intersect(DE_kegg, pw_cmpds)
    k <- length(overlap)
    n <- length(pw_cmpds)
    pval <- phyper(k - 1, n, bg_n - n, de_n, lower.tail = FALSE)
    data.frame(Pathway_ID = pid,
               Pathway_Name = pathway_list[pid],
               DE_in_Pathway = k,
               Pathway_Size = n,
               Pvalue = pval,
               Overlap = paste(overlap, collapse = ";"),
               stringsAsFactors = FALSE)
  })
  msea_df <- bind_rows(msea_list) %>%
    mutate(FDR = p.adjust(Pvalue, method = "BH"),
           RichFactor = DE_in_Pathway / Pathway_Size) %>%
    arrange(Pvalue)

  write.csv(msea_df, file.path(OUT_DIR, "metabolomics_KEGG_MSEA.csv"), row.names = FALSE)
  cat("  Saved: ", file.path(OUT_DIR, "metabolomics_KEGG_MSEA.csv"), "\n")

  # Bubble plot of top pathways (Pvalue < 0.05)
  top_msea <- msea_df %>%
    filter(Pvalue < 0.05) %>%
    mutate(Pathway_Name = str_remove(Pathway_Name, " - Homo sapiens \\(human\\)"),
           logP = -log10(Pvalue),
           Pathway_Name = factor(Pathway_Name, levels = rev(Pathway_Name)))

  if (nrow(top_msea) > 0) {
    p <- ggplot(top_msea, aes(x = RichFactor, y = Pathway_Name)) +
      geom_point(aes(size = DE_in_Pathway, color = logP)) +
      scale_color_gradient(low = "blue", high = "red", name = expression(-Log[10](P))) +
      scale_size_continuous(range = c(3, 10), name = "DE count") +
      labs(title = "Metabolite Set Enrichment (MSEA)",
           x = "Rich Factor", y = "") +
      theme_bw(base_size = 14) +
      theme(axis.text.y = element_text(size = 12, color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(OUT_DIR, "metabolomics_KEGG_bubble.pdf"), p, width = 10, height = 7)
    cat("  Saved: ", file.path(OUT_DIR, "metabolomics_KEGG_bubble.pdf"), "\n")
  } else {
    cat("  No pathways with P < 0.05.\n")
  }
} else {
  cat("  Too few KEGG IDs (< 3) for enrichment. Skipping MSEA.\n")
}

# ===== PART 5: GSEA (Gene-level, supplement) =====
cat("\n===== PART 5: GSEA (Gene Ontology) =====\n")
# Note: standard GSEA uses gene-level data. If RNA-seq expression data is
# available, load it here and run run_gsea(). Otherwise this section uses
# metabolite-level log2FC as a placeholder for gene set enrichment.
#
# Uncomment and adapt if a gene-level expression matrix is available:
#   rna_data <- read.csv(file.path(DATA_DIR, "RNA_expression.csv"), row.names = 1)
#   # ... run limma to get gene-level logFC, then:
#   geneList <- sort(setNames(gene_deg$logFC, gene_deg$ENTREZID), decreasing = TRUE)
#   gsea_res <- run_gsea(geneList)
#   save_enrichment(gsea_res$kegg, file.path(OUT_DIR, "metabolomics_GSEA_KEGG.csv"))
#   save_enrichment(gsea_res$go,   file.path(OUT_DIR, "metabolomics_GSEA_GO.csv"))

cat("  GSEA requires gene-level expression data (not loaded here).\n")
cat("  To enable: load RNA-seq counts, run limma DE, pass geneList to run_gsea().\n")

# ===== PART 6: CIBERSORT (Immune Infiltration) =====
cat("\n===== PART 6: CIBERSORT Immune Infiltration =====\n")
# CIBERSORT estimates immune cell fractions from bulk gene expression.
# Expected input: CIBERSORT output CSV with rows = samples, cols = cell types.
#
# Expected format:
#   cibersort_results.csv   (columns: Sample, B cells naive, ... , Neutrophils, P-value)
cibersort_file <- file.path(DATA_DIR, "cibersort_results.csv")

if (file.exists(cibersort_file)) {
  cib <- read.csv(cibersort_file, header = TRUE, row.names = 1)
  cib$Group <- ifelse(grepl("^TIGIT", rownames(cib)), "alphaTIGIT", "IgG")

  # Immune cell fraction heatmap
  frac_cols <- setdiff(colnames(cib), "Group")
  annotation_row <- data.frame(Group = cib$Group, row.names = rownames(cib))
  ann_cib_colors <- list(Group = c(IgG = "#4DBBD5", alphaTIGIT = "#E64B35"))

  pheatmap(t(as.matrix(cib[, frac_cols])),
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = annotation_row,
           annotation_colors = ann_cib_colors,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           cellwidth = 20, cellheight = 12,
           fontsize = 11,
           main = "CIBERSORT Immune Infiltration",
           filename = file.path(OUT_DIR, "cibersort_heatmap.pdf"))
  cat("  Saved: ", file.path(OUT_DIR, "cibersort_heatmap.pdf"), "\n")

  # Boxplot: compare each cell type between groups
  cib_long <- cib %>%
    rownames_to_column("Sample") %>%
    pivot_longer(all_of(frac_cols), names_to = "CellType", values_to = "Fraction")

  p <- ggplot(cib_long, aes(x = Group, y = Fraction, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c(IgG = "#4DBBD5", alphaTIGIT = "#E64B35")) +
    labs(title = "CIBERSORT: Immune Cell Fractions", x = NULL, y = "Fraction") +
    theme_classic() +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none") +
    geom_signif(comparisons = list(c("IgG", "alphaTIGIT")),
                test = "t.test", map_signif_level = FALSE,
                textsize = 3, tip_length = 0.02)
  ggsave(file.path(OUT_DIR, "cibersort_boxplots.pdf"), p,
         width = 14, height = ceiling(length(frac_cols) / 4) * 3.5)
  cat("  Saved: ", file.path(OUT_DIR, "cibersort_boxplots.pdf"), "\n")
} else {
  cat("  CIBERSORT file not found:", cibersort_file, "\n")
  cat("  Place cibersort_results.csv in", DATA_DIR, "to enable this section.\n")
}

# ===== SUMMARY =====
cat("\n========================================\n")
cat("  TIGIT Metabolomics Analysis Complete\n")
cat("========================================\n")
cat("Output files:\n")
cat("  ", file.path(OUT_DIR, "metabolomics_DE_results.csv"), "\n")
cat("  ", file.path(OUT_DIR, "metabolomics_PCA.pdf"), "\n")
cat("  ", file.path(OUT_DIR, "metabolomics_volcano.pdf"), "\n")
cat("  ", file.path(OUT_DIR, "metabolomics_heatmap_sig.pdf"), "\n")
cat("  ", file.path(OUT_DIR, "metabolomics_heatmap_top50.pdf"), "\n")
cat("  ", file.path(OUT_DIR, "metabolomics_boxplots.pdf"), "\n")
cat("  ", file.path(OUT_DIR, "metabolomics_KEGG_MSEA.csv"), "\n")
cat("  ", file.path(OUT_DIR, "metabolomics_KEGG_bubble.pdf"), "\n")
if (file.exists(cibersort_file)) {
  cat("  ", file.path(OUT_DIR, "cibersort_heatmap.pdf"), "\n")
  cat("  ", file.path(OUT_DIR, "cibersort_boxplots.pdf"), "\n")
}
cat("=== TIGIT Metabolomics Complete ===\n")
