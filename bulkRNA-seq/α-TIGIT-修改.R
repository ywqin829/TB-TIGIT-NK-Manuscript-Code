# =============================================================================
# alpha-TIGIT-modify.R -- Bulk RNA-seq: anti-TIGIT vs IgG in NK cells
# TB-TIGIT-NK manuscript
# Input:  DATA_DIR/raw_data.csv  |  Output: ./output/
# Panels: 1.Normalization+DE [if(FALSE)] 2.Volcano+Heatmap 3.GO/KEGG
#         4.GSEA  5.CIBERSORT  6.Spearman(TF vs phenotype)  7.Pathway heatmaps
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
source("../utils/go_kegg.R")
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
# PART 3: GO/KEGG enrichment (active)
# =============================================================================
cat("=== Part 3: GO/KEGG enrichment ===\n")
sig_csv <- file.path(OUT_DIR, "sig_gene.csv")
if (file.exists(sig_csv)) { sig_genes_symbol <- read.csv(sig_csv)$Gene
} else { sig_genes_symbol <- rownames(deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1.0, ]) }
cat("  Significant genes:", length(sig_genes_symbol), "\n")

if (length(sig_genes_symbol) >= 5) {
  # GO enrichment
  cat("  GO enrichment...\n")
  ego <- run_go_enrichment(sig_genes_symbol, ont = "all")
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    save_enrichment(ego, file.path(OUT_DIR, "ego.csv"))
    # Dotplot by ontology
    pdf(file.path(OUT_DIR, "go_dotplot.pdf"), width = 10, height = 8)
    print(dotplot(ego, split = "ONTOLOGY", showCategory = 15) + facet_grid(ONTOLOGY ~ ., scales = "free") +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) + theme(axis.text.y = element_text(size = 8)))
    dev.off()
    # Barplot
    pdf(file.path(OUT_DIR, "go_barplot.pdf"), width = 10, height = 8)
    print(barplot(ego, split = "ONTOLOGY", showCategory = 15) + facet_grid(ONTOLOGY ~ ., scales = "free") +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 60)))
    dev.off()
    # Network plots
    gl <- deg$logFC; names(gl) <- rownames(deg)
    pdf(file.path(OUT_DIR, "go_cnet.pdf"), width = 12, height = 10)
    print(cnetplot(ego, foldChange = gl, showCategory = 6, layout = "gem")); dev.off()
    pdf(file.path(OUT_DIR, "go_emap.pdf"), width = 10, height = 8)
    print(emapplot(enrichplot::pairwise_termsim(ego), showCategory = 15, layout = "kk", cex_category = 1.5, min_edge = 0.8)); dev.off()
    pdf(file.path(OUT_DIR, "go_heatplot.pdf"), width = 12, height = 8)
    print(heatplot(ego, foldChange = gl) + theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 4))); dev.off()
  }
  # KEGG enrichment
  cat("  KEGG enrichment...\n")
  options(timeout = 600)
  kk <- run_kegg_enrichment(sig_genes_symbol)
  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    save_enrichment(kk, file.path(OUT_DIR, "kegg.csv"))
    pdf(file.path(OUT_DIR, "kegg_dotplot.pdf"), width = 9, height = 7)
    print(dotplot(kk, showCategory = 20) + scale_y_discrete(labels = function(x) str_wrap(x, width = 50)))
    dev.off()
    # RichFactor dotplot
    kk_df <- as.data.frame(kk)
    if (nrow(kk_df) > 0) {
      kk_f <- kk_df %>% arrange(desc(RichFactor)) %>% mutate(Description = factor(Description, levels = rev(unique(Description))))
      p <- ggplot(kk_f, aes(x = RichFactor, y = Description, color = p.adjust, size = Count)) + geom_point() +
        scale_color_gradient(low = "red", high = "blue") + scale_size_continuous(range = c(3, 10)) +
        labs(title = "KEGG by RichFactor", x = "RichFactor", y = "") + theme_bw(base_size = 14) +
        theme(axis.text.y = element_text(size = 12, colour = "black"), plot.title = element_text(hjust = 0))
      ggsave(file.path(OUT_DIR, "kegg_richfactor.pdf"), p, width = 9, height = 7)
    }
  }
} else { cat("  WARNING: Too few genes (< 5) for enrichment\n") }
cat("=== Part 3 complete ===\n\n")

# =============================================================================
# PART 4: GSEA (active)
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

# =============================================================================
# PART 5: CIBERSORT immune infiltration (active)
# =============================================================================
cat("=== Part 5: CIBERSORT ===\n")
if (!requireNamespace("CIBERSORT", quietly = TRUE)) {
  cat("  NOTE: CIBERSORT not installed. Install via: devtools::install_github('Moonerss/CIBERSORT')\n")
} else {
  library(CIBERSORT); data(LM22)
  if (exists("exp_vst")) { mm <- exp_vst } else { mm <- log2(exp + 1) }
  cat("  Running CIBERSORT...\n")
  ci_res <- tryCatch(cibersort(sig_matrix = LM22, mixture_file = as.data.frame(mm), perm = 100, QN = FALSE),
    error = function(e) { cat("  CIBERSORT error:", conditionMessage(e), "\n"); return(NULL) })
  if (!is.null(ci_res)) {
    ci_all <- tibble::rownames_to_column(as.data.frame(ci_res), var = "Sample")
    write.csv(ci_all, file.path(OUT_DIR, "cibersort_all.csv"), row.names = FALSE)
    ci_dat <- ci_all[, 1:23]; colnames(ci_dat)[1] <- "Sample"
    # 5a. Boxplot
    ci_long <- ci_dat %>% left_join(group_info, by = "Sample") %>%
      reshape2::melt(id.vars = c("Sample", "Group"), variable.name = "celltype", value.name = "composition")
    p_box <- ggplot(ci_long, aes(x = celltype, y = composition)) + geom_boxplot(aes(fill = Group), position = position_dodge(0.5), width = 0.5) +
      labs(y = "Cell composition", x = "") + scale_fill_manual(values = c("IgG" = "#1f77b4", "alphaTIGIT" = "#ff7f0e")) + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), legend.position = "top") +
      stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test", hide.ns = TRUE)
    ggsave(file.path(OUT_DIR, "cibersort_boxplot.pdf"), p_box, height = 7, width = 10)
    # 5b. Stacked barplot
    ci_plot <- reshape2::melt(ci_dat, id.vars = "Sample", variable.name = "celltype", value.name = "composition")
    cols <- c(RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(12, "Set3"))
    p_bar <- ggplot(ci_plot, aes(x = Sample, y = composition, fill = celltype)) + geom_bar(position = "stack", stat = "identity") +
      scale_fill_manual(values = cols) + labs(x = "", y = "Estimated Proportion", title = "CIBERSORT") +
      scale_y_continuous(expand = c(0, 0)) + guides(fill = guide_legend(ncol = 1)) + theme_bw() +
      theme(legend.key = element_blank(), legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUT_DIR, "cibersort_barplot.pdf"), p_bar, height = 8, width = 8)
    # 5c. Heatmap
    k <- apply(ci_res, 2, function(x) sum(x == 0) < nrow(ci_res) / 2)
    if (sum(k) >= 2) {
      re2 <- as.data.frame(t(ci_res[, k]))
      an <- data.frame(Group = group_info[colnames(exp), "Group"], row.names = colnames(exp))[colnames(re2), , drop = FALSE]
      pdf(file.path(OUT_DIR, "cibersort_heatmap.pdf"), width = 7, height = 6)
      pheatmap(re2, scale = "row", show_colnames = FALSE, cluster_cols = FALSE, annotation_col = an,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100), breaks = seq(-3, 3, length.out = 100),
        main = "CIBERSORT Cell Fractions")
      dev.off()
    }
  }
}
cat("=== Part 5 complete ===\n\n")

# =============================================================================
# PART 6: Spearman correlation (TF vs phenotype) (active)
# =============================================================================
cat("=== Part 6: Spearman (TF vs phenotype) ===\n")
if (file.exists(exp_vst_csv)) {
  if (!exists("exp_vst")) exp_vst <- as.matrix(read.csv(exp_vst_csv, row.names = 1))
  all_genes <- c("ZEB2", "LITAF", "EOMES", "TBX21", "PRDM1", "ZBED2", "FOS", "ESR2", "SMAD9",
                 "TNF", "GZMA", "GZMB", "MPEG1", "MYO10", "RAC1", "CDC42", "RHOA")
  avail <- intersect(all_genes, rownames(exp_vst))
  cat("  Available:", length(avail), "/", length(all_genes), "\n")
  if (length(avail) >= 4) {
    cr <- corr.test(t(exp_vst[avail, , drop = FALSE]) %>% as.data.frame(), method = "spearman", adjust = "fdr")
    pdf(file.path(OUT_DIR, "spearman_tf_phenotype.pdf"), width = 8, height = 7)
    corrplot(cr$r, method = "pie", type = "lower", tl.col = "black", tl.srt = 45, tl.cex = 0.7,
      p.mat = cr$p, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", pch.cex = 0.8, pch.col = "white",
      col = colorRampPalette(c("blue", "white", "red"))(50), title = "TF vs Phenotype (Spearman)", mar = c(1, 1, 2, 1))
    dev.off()
    write.csv(cr$r, file.path(OUT_DIR, "spearman_cor_matrix.csv"))
  } else { cat("  WARNING: Too few genes for correlation\n") }
} else { cat("  WARNING: exp_vst.csv missing, skipping Part 6\n") }
cat("=== Part 6 complete ===\n\n")

# =============================================================================
# PART 7: Key metabolic pathway heatmaps (active)
# =============================================================================
cat("=== Part 7: Metabolic pathway heatmaps ===\n")
# Helper functions
extract_pathway_expr <- function(expr, genes, name) {
  avail <- intersect(genes, rownames(expr))
  if (length(avail) == 0) { cat("  --", name, ": no genes\n"); return(NULL) }
  cat("  --", name, ":", length(avail), "/", length(genes), "\n"); expr[avail, , drop = FALSE]
}
plot_pathway_heatmap <- function(mat, name, grp) {
  ms <- t(scale(t(as.matrix(mat))))
  cc <- intersect(colnames(ms), rownames(grp))
  ms <- ms[, cc, drop = FALSE]
  ht <- Heatmap(ms, name = "Z-Score", col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
    top_annotation = HeatmapAnnotation(Group = grp[cc, "Group"],
      col = list(Group = c("IgG" = "#1f77b4", "alphaTIGIT" = "#ff7f0e"))),
    cluster_columns = FALSE, column_split = factor(grp[cc, "Group"], levels = c("IgG", "alphaTIGIT")),
    column_gap = unit(3, "mm"), cluster_rows = TRUE, show_row_names = TRUE, show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 10), column_title = name, column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  draw(ht)
}
# Pathway gene lists
pathways <- list(
  `HIF-1a_Pathway` = c("HIF1A", "SLC2A1", "LDHA", "PDK1", "ENO1", "PFKFB3", "PKM", "HK2", "ALDOA", "PGK1"),
  `Polyamine_Pathway` = c("ODC1", "AMD1", "SAT1", "SMS", "SMOX", "ARG1", "ARG2", "OAZ1", "AZIN1"),
  `Pyrimidine_Pathway` = c("CAD", "DHODH", "UMPS", "UCK2", "TYMS", "TK1", "DTYMK", "RRM1", "RRM2"),
  `Glutamine_Pathway` = c("GLS", "GLS2", "GLUD1", "GLUL", "SLC1A5", "SLC7A5", "SLC3A2", "GOT1", "GOT2"),
  `mTOR_Pathway` = c("MTOR", "RPTOR", "RICTOR", "RPS6KB1", "EIF4EBP1", "AKT1", "AKT2", "TSC1", "TSC2", "RHEB"),
  `NK_Effector` = c("GZMB", "PRF1", "IFNG", "TNF", "FASLG", "GNLY", "GZMA", "GZMH", "GZMK", "GZMM"),
  `Epigenetic_Enzymes` = c("HDAC1", "HDAC2", "HDAC3", "HDAC6", "KAT2A", "KAT2B", "EP300", "CREBBP", "EZH2", "KDM5A", "KDM6A")
)
# Expression matrix
if (file.exists(exp_vst_csv) && !exists("exp_vst")) exp_vst <- as.matrix(read.csv(exp_vst_csv, row.names = 1))
expr_mat <- if (exists("exp_vst")) exp_vst else log2(exp + 1)
cat("  Expression matrix:", nrow(expr_mat), "genes\n")
# Generate heatmaps
dir.create(file.path(OUT_DIR, "pathway_heatmaps"), showWarnings = FALSE)
for (pn in names(pathways)) {
  e <- extract_pathway_expr(expr_mat, pathways[[pn]], pn)
  if (!is.null(e)) {
    pdf(file.path(OUT_DIR, "pathway_heatmaps", paste0(pn, "_heatmap.pdf")), width = 8, height = max(4, nrow(e) * 0.4 + 2))
    plot_pathway_heatmap(e, gsub("_", " ", pn), group_info); dev.off()
  }
}
cat("=== Part 7 complete ===\n\n")

# ===== SUMMARY =====
cat("========================================\n=== ANALYSIS COMPLETE ===\n========================================\n")
cat("Output:", normalizePath(OUT_DIR, mustWork = FALSE), "\n\nFiles:\n")
for (f in sort(list.files(OUT_DIR, recursive = TRUE))) cat("  ", file.path(OUT_DIR, f), "\n")
cat("\n=== End of script ===\n")
