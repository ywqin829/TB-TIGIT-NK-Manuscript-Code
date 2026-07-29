# =============================================================================
# alpha-TIGIT-modify.R -- Bulk RNA-seq: anti-TIGIT vs IgG in NK cells
# TB-TIGIT-NK manuscript
# Input:  DATA_DIR/raw_data.csv  |  Output: ./output/
# Panels: 1.Normalization+DE [if(FALSE)] 2.Volcano 3.GSEA
# =============================================================================

# ===== SETUP =====
rm(list = ls()); cat("=== alpha-TIGIT Bulk RNA-seq Analysis ===\n")
suppressPackageStartupMessages({
  library(limma); library(edgeR); library(DESeq2); library(dplyr)
  library(tidyr); library(ggplot2); library(pheatmap); library(ggrepel)
  library(clusterProfiler); library(org.Hs.eg.db); library(EnhancedVolcano)
  library(factoextra); library(GSVA); library(tibble); library(stringr)
  library(corrplot); library(psych); library(ComplexHeatmap); library(circlize)
  library(reshape2); library(ggpubr)
})
source("utils/go_kegg.R")
OUT_DIR <- "./output"; dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
DATA_DIR <- "./data"
cat("Output:", normalizePath(OUT_DIR, mustWork = FALSE),
    "  Data:", normalizePath(DATA_DIR, mustWork = FALSE), "\n\n")

# ===== DATA LOADING =====
cat("=== Loading raw data ===\n")
raw_csv <- file.path(DATA_DIR, "raw_data.csv")
if (!file.exists(raw_csv)) stop("raw_data.csv not found at: ", normalizePath(raw_csv, mustWork = FALSE))
exp_raw <- read.csv(raw_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
if ("Symbol" %in% colnames(exp_raw)) { rownames(exp_raw) <- make.unique(exp_raw$Symbol); exp_raw$Symbol <- NULL
} else if ("symbol" %in% colnames(exp_raw)) { rownames(exp_raw) <- make.unique(exp_raw$symbol); exp_raw$symbol <- NULL }
cat("  Raw dimensions:", nrow(exp_raw), "genes x", ncol(exp_raw), "samples\n")
# Filter low-expressed genes and "-." artefacts
exp <- exp_raw[rowSums(exp_raw != 0) > 1, ]
exp <- exp[!grepl("^-|\\.-", rownames(exp)), ]
cat("  After filtering:", nrow(exp), "genes x", ncol(exp), "samples\n")
# Build group info from sample names
samples <- colnames(exp)
group_info <- data.frame(Sample = samples,
  Group = factor(ifelse(grepl("^Ctrl", samples), "IgG",
                 ifelse(grepl("^TIGIT", samples), "alphaTIGIT", NA)),
                 levels = c("IgG", "alphaTIGIT")), row.names = samples)
group_info <- subset(group_info, Group %in% c("IgG", "alphaTIGIT") & !is.na(Group))
common_samples <- intersect(colnames(exp), rownames(group_info))
exp <- exp[, common_samples]; group_info <- group_info[common_samples, , drop = FALSE]
group_list <- group_info$Group; names(group_list) <- rownames(group_info)
cat("  Groups:\n"); print(table(group_list)); cat("\n")

# =============================================================================
# PART 1: Normalization + limma DE [if (FALSE) -- run once]
# =============================================================================
if (FALSE) {
  cat("=== Part 1: Normalization + limma DE ===\n")
  # 1a. VST normalization
  cat("  VST normalization...\n")
  dds <- DESeqDataSetFromMatrix(countData = round(exp), colData = group_info, design = ~ Group)
  dds <- estimateSizeFactors(dds)
  exp_vst <- assay(vst(dds, blind = TRUE, fitType = "local"))
  # 1b. Boxplots
  colors_box <- c("IgG" = "green", "alphaTIGIT" = "blue")
  pdf(file.path(OUT_DIR, "boxplot_raw.pdf"), width = 8, height = 5)
  par(mar = c(7, 4, 2, 2))
  boxplot(exp, outline = FALSE, notch = FALSE, col = colors_box[as.character(group_list)], las = 2, main = "Raw counts", ylab = "Count")
  dev.off()
  pdf(file.path(OUT_DIR, "boxplot_vst.pdf"), width = 8, height = 5)
  par(mar = c(7, 4, 2, 2))
  boxplot(exp_vst, outline = FALSE, notch = FALSE, col = colors_box[as.character(group_list)], las = 2, main = "Expression (VST)", ylab = "VST")
  dev.off()
  # 1c. PCA
  cat("  PCA...\n")
  dat_pca <- t(exp_vst); dat_pca <- dat_pca[, apply(dat_pca, 2, sd) != 0]
  pca_res <- prcomp(dat_pca, scale. = TRUE)
  expl_var <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
  p_pca <- fviz_pca_ind(pca_res, geom.ind = "point", col.ind = group_list,
    palette = c("IgG" = "green", "alphaTIGIT" = "blue"),
    addEllipses = TRUE, ellipse.type = "confidence", legend.title = "Group") +
    labs(title = "PCA", x = paste0("PC1 (", round(expl_var[1], 1), "%)"), y = paste0("PC2 (", round(expl_var[2], 1), "%)"))
  ggsave(file.path(OUT_DIR, "pca.pdf"), p_pca, width = 6, height = 5)
  ggsave(file.path(OUT_DIR, "pca_scree.pdf"), fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 50)), width = 5, height = 4)
  write.csv(exp_vst, file.path(OUT_DIR, "exp_vst.csv"), row.names = TRUE)
  # 1d. limma-voom DE
  cat("  limma-voom DE...\n")
  deg_obj <- DGEList(counts = exp, group = group_list)
  deg_obj <- deg_obj[filterByExpr(deg_obj, group = group_list), , keep.lib.sizes = FALSE]
  deg_obj <- calcNormFactors(deg_obj)
  design <- model.matrix(~ 0 + group_list); colnames(design) <- levels(group_list)
  v <- voom(deg_obj, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, makeContrasts(alphaTIGIT_vs_IgG = alphaTIGIT - IgG, levels = design))
  fit2 <- eBayes(fit2)
  deg <- topTable(fit2, coef = "alphaTIGIT_vs_IgG", number = Inf, adjust.method = "BH", sort.by = "P")
  logFC_cut <- 1.0; pval_cut <- 0.05
  deg$change <- ifelse(deg$adj.P.Val < pval_cut & deg$logFC < -logFC_cut, "down",
                ifelse(deg$adj.P.Val < pval_cut & deg$logFC > logFC_cut, "up", "stable"))
  cat("  DEGs:", nrow(deg), "| up:", sum(deg$change == "up"), "down:", sum(deg$change == "down"), "\n")
  write.csv(cbind(Gene = rownames(deg), deg), file.path(OUT_DIR, "deg_limma_voom.csv"), row.names = FALSE)
  write.csv(deg, file.path(OUT_DIR, "deg.csv"), row.names = TRUE)
  up_genes <- rownames(deg[deg$change == "up", ]); down_genes <- rownames(deg[deg$change == "down", ])
  write.table(up_genes, file.path(OUT_DIR, "up_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(down_genes, file.path(OUT_DIR, "down_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  sig_gene <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1.0, ]
  write.csv(cbind(Gene = rownames(sig_gene), sig_gene), file.path(OUT_DIR, "sig_gene.csv"), row.names = FALSE)
  cat("=== Part 1 complete ===\n\n")
}

# =============================================================================
# PART 2: Volcano plot + DEG heatmap (active)
# =============================================================================
cat("=== Part 2: Volcano plot + DEG heatmap ===\n")
deg_csv <- file.path(OUT_DIR, "deg.csv")
if (!file.exists(deg_csv)) stop("Run Part 1 first to generate: ", deg_csv)
deg <- read.csv(deg_csv, row.names = 1)
cat("  DEGs loaded:", nrow(deg), "\n")
if (!"change" %in% colnames(deg)) {
  deg$change <- ifelse(deg$adj.P.Val < 0.05 & deg$logFC < -1, "down",
                ifelse(deg$adj.P.Val < 0.05 & deg$logFC > 1, "up", "stable"))
}
# Volcano
p_volc <- ggplot(deg, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.6, size = 2.0, aes(color = change)) +
  scale_color_manual(values = c("down" = "#0000EE", "stable" = "grey", "up" = "#8B2323")) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
  ylab("-Log10(P.Value)") + theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 15), legend.text = element_text(size = 14))
gene_show <- c("PLIN5", "CPT1A", "PPARG", "APOC2", "MLXIPL", "INSR", "IRS2", "RIPK3", "GPX3", "MAPT", "COX5B1", "ZBTB20", "IL1B")
deg$Gene <- rownames(deg)
label_df <- deg[deg$Gene %in% gene_show, ]
p_volcano <- p_volc + geom_point(size = 3, shape = 1, data = label_df) +
  geom_label_repel(data = label_df, aes(label = Gene), color = "black", size = 5, max.overlaps = Inf)
ggsave(file.path(OUT_DIR, "volcano.pdf"), p_volcano, width = 8, height = 7)
ggsave(file.path(OUT_DIR, "volcano.png"), p_volcano, width = 8, height = 7, dpi = 300)
# DEG heatmap
exp_vst_csv <- file.path(OUT_DIR, "exp_vst.csv")
if (file.exists(exp_vst_csv)) {
  exp_vst <- as.matrix(read.csv(exp_vst_csv, row.names = 1))
  deg_genes <- intersect(rownames(deg[deg$change != "stable", ]), rownames(exp_vst))
  if (length(deg_genes) >= 2) {
    ann_col <- data.frame(Group = group_info[colnames(exp_vst), "Group"], row.names = colnames(exp_vst))
    pdf(file.path(OUT_DIR, "deg_heatmap.pdf"), width = 8, height = 7)
    pheatmap(exp_vst[deg_genes, ], annotation_col = ann_col,
      annotation_colors = list(Group = c("IgG" = "#1f77b4", "alphaTIGIT" = "#ff7f0e")),
      scale = "row", show_rownames = FALSE, show_colnames = FALSE,
      color = colorRampPalette(c("navy", "cyan", "white", "red", "firebrick"))(400), main = paste0("DEGs (", length(deg_genes), ")"))
    dev.off()
  }
}
cat("=== Part 2 complete ===\n\n")

# =============================================================================
# PART 3: GSEA (active)
# =============================================================================
cat("=== Part 4: GSEA ===\n")
deg$Gene <- rownames(deg)
gene_map <- tryCatch(
  AnnotationDbi::select(org.Hs.eg.db, keys = as.character(deg$Gene), keytype = "SYMBOL", columns = c("ENTREZID", "SYMBOL")),
  error = function(e) { cat("  WARNING: gene mapping failed:", conditionMessage(e), "\n"); return(NULL) })
if (!is.null(gene_map) && nrow(gene_map) > 0) {
  gene_map <- gene_map %>% distinct(SYMBOL, .keep_all = TRUE) %>% filter(!is.na(ENTREZID))
  colnames(gene_map)[colnames(gene_map) == "SYMBOL"] <- "Gene"
  temp <- inner_join(gene_map, deg, by = "Gene") %>% dplyr::select(ENTREZID, logFC) %>% distinct(ENTREZID, .keep_all = TRUE)
  geneList <- temp$logFC; names(geneList) <- as.character(temp$ENTREZID)
  geneList <- geneList[!is.na(names(geneList))]
  set.seed(42); geneList <- geneList + runif(length(geneList), min = -1e-5, max = 1e-5); geneList <- sort(geneList, decreasing = TRUE)
  cat("  Gene list size:", length(geneList), "\n")
  gsea_res <- run_gsea(geneList)
  # GO GSEA
  if (!is.null(gsea_res$go) && nrow(as.data.frame(gsea_res$go)) > 0) {
    save_enrichment(gsea_res$go, file.path(OUT_DIR, "gsea_go.csv"))
    top_ids <- gsea_res$go@result$ID[1:min(10, nrow(gsea_res$go@result))]
    dir.create(file.path(OUT_DIR, "gsea_go_plots"), showWarnings = FALSE)
    for (gid in top_ids) {
      desc <- gsea_res$go@result$Description[gsea_res$go@result$ID == gid]
      if (length(desc) == 1) { pdf(file.path(OUT_DIR, "gsea_go_plots", paste0(gsub("/", "_", gid), ".pdf")), width = 7, height = 5)
        print(gseaplot2(gsea_res$go, geneSetID = gid, title = desc, pvalue_table = TRUE)); dev.off() }
    }
  }
  # KEGG GSEA
  if (!is.null(gsea_res$kegg) && nrow(as.data.frame(gsea_res$kegg)) > 0) {
    save_enrichment(gsea_res$kegg, file.path(OUT_DIR, "gsea_kegg.csv"))
    top_ids <- gsea_res$kegg@result$ID[1:min(5, nrow(gsea_res$kegg@result))]
    dir.create(file.path(OUT_DIR, "gsea_kegg_plots"), showWarnings = FALSE)
    for (kid in top_ids) {
      desc <- gsea_res$kegg@result$Description[gsea_res$kegg@result$ID == kid]
      if (length(desc) == 1) { pdf(file.path(OUT_DIR, "gsea_kegg_plots", paste0(kid, ".pdf")), width = 7, height = 5)
        print(gseaplot2(gsea_res$kegg, geneSetID = kid, title = desc, pvalue_table = TRUE)); dev.off() }
    }
  }
  write.csv(data.frame(ENTREZID = names(geneList), logFC = geneList), file.path(OUT_DIR, "gsea_geneList.csv"), row.names = FALSE)
} else { cat("  WARNING: Gene mapping failed, skipping GSEA\n") }
cat("=== Part 4 complete ===\n\n")

# ===== SUMMARY =====
cat("========================================\n=== ANALYSIS COMPLETE ===\n========================================\n")
cat("Output:", normalizePath(OUT_DIR, mustWork = FALSE), "\n\nFiles:\n")
for (f in sort(list.files(OUT_DIR, recursive = TRUE))) cat("  ", file.path(OUT_DIR, f), "\n")
cat("\n=== End of script ===\n")
