# ==============================================================================
# TIGIT代谢组学测序分析.R — Complete Metabolomics + Multi-omics Integration
# TB-TIGIT-NK manuscript
#
# Input:
#   ./data/111.csv              — Metabolomics expression matrix (metabolites x samples)
#   ./data/name_map_1.csv       — Metabolite name → KEGG/HMDB ID mapping
#   ./data/msea_ora.csv         — (optional) MetaboAnalyst ORA results
#   ../bulkRNA-seq/data/exp_vst — RNA-seq VST expression (for DIABLO / Part 10)
#   ../bulkRNA-seq/data/sig_gene.csv — DEG table (for Pathview / Part 8)
#
# Output (Fig 6 panels):
#   Fig 6A — Volcano plot (28 DEMs, VIP > 1, P < 0.05)
#   Fig 6B — KEGG pathway enrichment bubble (MSEA)
#   Fig 6D — Boxplots of 7 key metabolites + Pathview KEGG pathway maps
#   Fig 6E — DIABLO multi-omics integration (mixOmics block.splsda)
# ==============================================================================

# ===== SETUP =====
rm(list = ls())
cat("=== TIGIT Metabolomics + Multi-omics Integration ===\n")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(ggrepel); library(ggsignif); library(pheatmap); library(reshape2)
  library(FactoMineR); library(factoextra); library(ropls)
  library(KEGGREST); library(clusterProfiler); library(enrichplot)
  library(pathview); library(org.Hs.eg.db)
  library(psych); library(corrplot); library(randomForest)
  library(mixOmics); library(stringr)
})

OUT_DIR <- "./output"
DATA_DIR <- "./data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# ===== LOAD DATA =====
cat("\n===== Loading metabolomics data =====\n")
metadata <- read.csv(file.path(DATA_DIR, "111.csv"), header = TRUE, row.names = 1)
cat(sprintf("  Raw: %d metabolites x %d samples\n", nrow(metadata), ncol(metadata)))

meta_exp <- metadata[rowSums(metadata != 0) > 1, ]
cat(sprintf("  After filtering: %d metabolites\n", nrow(meta_exp)))

samples <- colnames(meta_exp)
group_info <- data.frame(
  Sample = samples,
  Group = ifelse(grepl("^IgG", samples), "IgG",
         ifelse(grepl("^TIGIT", samples), "alphaTIGIT", NA))
)
group_info <- subset(group_info, Group %in% c("IgG", "alphaTIGIT"))
group_info <- group_info[order(factor(group_info$Group, levels = c("IgG", "alphaTIGIT"))), ]
group_list <- factor(group_info$Group, levels = c("IgG", "alphaTIGIT"))
cat(sprintf("  IgG: %d, alphaTIGIT: %d\n", sum(group_list == "IgG"), sum(group_list == "alphaTIGIT")))

# =============================================================================
# PART 1: PCA
# =============================================================================
cat("\n===== PART 1: PCA =====\n")
pca_res <- prcomp(t(meta_exp), scale. = TRUE)
explained_var <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
p_pca <- fviz_pca_ind(pca_res, geom.ind = "point", col.ind = group_info$Group,
             palette = c("green", "blue"), addEllipses = TRUE,
             ellipse.type = "confidence", legend.title = "Group") +
  labs(title = "PCA (2D) of Metabolomics",
       x = paste0("PC1 (", round(explained_var[1], 1), "%)"),
       y = paste0("PC2 (", round(explained_var[2], 1), "%)"))
ggsave(file.path(OUT_DIR, "metabolomics_PCA.pdf"), p_pca, width = 7, height = 6)

# =============================================================================
# PART 2: Differential Expression + VIP (Fig 6A prep)
# =============================================================================
cat("\n===== PART 2: Differential Expression =====\n")

# FC, log2FC, P-value
res_df <- apply(meta_exp, 1, function(x) {
  x1 <- x[group_info$Sample[group_info$Group == "alphaTIGIT"]]
  x2 <- x[group_info$Sample[group_info$Group == "IgG"]]
  fc <- mean(x1) / mean(x2)
  log2fc <- log2(fc)
  pval <- t.test(x1, x2)$p.value
  return(c(FC = fc, log2FC = log2fc, P.Value = pval))
})
res_df <- as.data.frame(t(res_df))
res_df$Metabolite <- rownames(res_df)
res_df$FDR <- p.adjust(res_df$P.Value, method = "fdr")

# PLS-DA VIP scores
meta_exp_t <- t(meta_exp)
pls_model <- opls(meta_exp_t, group_info$Group, predI = 1, orthoI = 1, crossvalI = 6)
vip_scores <- getVipVn(pls_model)
res_df$VIP <- vip_scores[rownames(res_df)]

# Classification: VIP > 1 & P.Value < 0.05
res_df$type <- "insig"
res_df$type[res_df$VIP > 1 & res_df$P.Value < 0.05 & res_df$log2FC > 0] <- "up"
res_df$type[res_df$VIP > 1 & res_df$P.Value < 0.05 & res_df$log2FC < 0] <- "down"

res_df <- res_df[, c("Metabolite", "VIP", "P.Value", "FDR", "FC", "log2FC", "type")]
write.csv(res_df, file.path(OUT_DIR, "metabolite_statistics.csv"), row.names = FALSE)
cat(sprintf("  Up: %d, Down: %d\n", sum(res_df$type == "up"), sum(res_df$type == "down")))

# =============================================================================
# PART 3: Volcano Plot (Fig 6A)
# =============================================================================
cat("\n===== PART 3: Volcano Plot (Fig 6A) =====\n")

res_df$negLog10P <- -log10(res_df$P.Value)

sig_up   <- res_df %>% filter(type == "up")   %>% arrange(desc(log2FC), desc(negLog10P)) %>% head(20)
sig_down <- res_df %>% filter(type == "down") %>% arrange(log2FC, desc(negLog10P)) %>% head(20)
label_df <- rbind(sig_up, sig_down)

p_volcano <- ggplot(res_df, aes(x = log2FC, y = negLog10P, color = type, size = VIP)) +
  geom_point(alpha = 0.7, shape = 16) +
  scale_color_manual(values = c("up" = "#8B2323", "down" = "#0000EE", "insig" = "grey")) +
  scale_size_continuous(range = c(0.1, 6), breaks = c(0.5, 1, 1.5), name = "VIP") +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
  geom_point(data = label_df, aes(x = log2FC, y = negLog10P),
             shape = 1, size = 3, color = "black", inherit.aes = FALSE, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = label_df,
    aes(x = log2FC, y = negLog10P, label = Metabolite),
    size = 5, max.overlaps = 20, inherit.aes = FALSE) +
  theme_bw() +
  labs(title = "Volcano Plot: Significant Metabolites",
       subtitle = "VIP > 1 & P.Value < 0.05",
       x = "log2(Fold Change)", y = "-log10(P-value)") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14), axis.text = element_text(size = 12))

ggsave(file.path(OUT_DIR, "metabolomics_volcano.pdf"), p_volcano, width = 10, height = 8)

# =============================================================================
# PART 4: Boxplots of 7 Key Metabolites (Fig 6D)
# =============================================================================
cat("\n===== PART 4: Boxplots — 7 Key Metabolites (Fig 6D) =====\n")

sig_metabolites <- res_df$Metabolite[res_df$type != "insig"]
exp_sig <- meta_exp[sig_metabolites, ]

metabolite_names <- c('Beta-Alanine', '1,3-diaminopropane', 'Uracil',
                      'Spermidine', 'Glutamine', "5'-Cytidylic acid",
                      'Dephosphocoenzyme A (Dephospho-CoA)')
exp_sig_select <- exp_sig[intersect(metabolite_names, rownames(exp_sig)), ]
write.csv(exp_sig_select, file.path(OUT_DIR, "exp_sig_select.csv"), row.names = TRUE)

df_long <- exp_sig_select %>%
  as.data.frame() %>% rownames_to_column("Metabolite") %>%
  pivot_longer(-Metabolite, names_to = "Sample", values_to = "Concentration") %>%
  left_join(group_info, by = c("Sample"))

p_box <- ggplot(df_long, aes(x = Group, y = Concentration, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(), width = 0.6, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.6)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.6) +
  facet_wrap(~ Metabolite, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Concentration", title = "Concentration of Significant Metabolites") +
  scale_fill_discrete(guide = "none") + theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 15),
        axis.title = element_text(size = 20, face = 'bold'),
        strip.text = element_text(size = 16), legend.position = "none") +
  geom_signif(comparisons = list(c("IgG", "alphaTIGIT")),
              test = "t.test", map_signif_level = FALSE, textsize = 4)

ggsave(file.path(OUT_DIR, "metabolomics_7key_boxplots.pdf"), p_box, width = 16, height = 6)

# =============================================================================
# PART 5: Heatmaps
# =============================================================================
cat("\n===== PART 5: Heatmaps =====\n")

annotation_col <- data.frame(Group = group_list, row.names = group_info$Sample)

# Top 50 by |log2FC|
top_fc <- res_df %>% arrange(desc(abs(log2FC))) %>% slice(1:min(50, n())) %>% pull(Metabolite)
meta_exp_top_fc <- meta_exp[top_fc, group_info$Sample, drop = FALSE]
pdf(file.path(OUT_DIR, "metabolomics_heatmap_top50.pdf"), width = 10, height = 14)
pheatmap(meta_exp_top_fc, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_col = annotation_col,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 15, cellheight = 10, main = "Top 50 Variable Metabolites")
dev.off()

# Significant metabolites
meta_exp_sig <- meta_exp[sig_metabolites, group_info$Sample, drop = FALSE]
if (nrow(meta_exp_sig) > 1) {
  pdf(file.path(OUT_DIR, "metabolomics_heatmap_sig.pdf"), width = 10, height = max(5, nrow(meta_exp_sig) * 0.3))
  pheatmap(meta_exp_sig, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
           show_colnames = FALSE, annotation_col = annotation_col,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           cellwidth = 20, cellheight = 10, main = "Significant Metabolites (VIP>1, P<0.05)")
  dev.off()
}

# =============================================================================
# PART 6: KEGG Enrichment — MSEA (Fig 6B)
# =============================================================================
cat("\n===== PART 6: KEGG Enrichment (MSEA) =====\n")

# Map metabolite names to KEGG IDs
name_map <- read.csv(file.path(DATA_DIR, "name_map_1.csv"), header = TRUE, row.names = 1)
name_map <- name_map %>% as.data.frame() %>% tibble::rownames_to_column("Metabolite")
sig_df <- res_df[res_df$type != "insig", ]
sig_df_KEGG <- sig_df %>% left_join(name_map %>% dplyr::select(Metabolite, KEGG), by = "Metabolite")

get_kegg_id <- function(metabolite_name) {
  tryCatch({
    kegg_info <- keggFind("compound", metabolite_name)
    if (length(kegg_info) > 0 && !is.null(names(kegg_info)))
      return(sub("^cpd:", "", names(kegg_info)[1]))
    return(NA)
  }, error = function(e) return(NA))
}

# Use MetaboAnalystR ORA if pre-computed; otherwise run KEGGREST MSEA
msea_ora_path <- file.path(DATA_DIR, "msea_ora.csv")
if (file.exists(msea_ora_path)) {
  msea_ora <- read.csv(msea_ora_path, header = TRUE, row.names = 1)
  msea_ora$Pathway <- rownames(msea_ora)
  msea_ora <- msea_ora %>% mutate(EnrichmentRatio = hits / expected)
  msea_ora_filtered <- msea_ora %>%
    filter(Raw.p < 0.25) %>%
    mutate(Pathway = factor(Pathway, levels = rev(Pathway)), logP = -log10(Raw.p))

  p_enrich <- ggplot(msea_ora_filtered, aes(x = logP, y = Pathway, size = EnrichmentRatio, color = Raw.p)) +
    geom_point(alpha = 0.85) +
    geom_text(data = msea_ora_filtered %>% filter(FDR < 0.05),
              aes(label = sprintf("%.3f", FDR)), nudge_x = 0.1, size = 3.5, color = "black") +
    scale_color_gradient(low = "red", high = "blue", name = "P-value") +
    scale_size_continuous(range = c(3, 13), name = "Enrichment Ratio") +
    labs(title = "Metabolite Pathway Enrichment Analysis",
         x = expression(-Log[10](P-value)), y = "") +
    theme_light(base_size = 14)
  ggsave(file.path(OUT_DIR, "metabolomics_KEGG_bubble.pdf"), p_enrich, width = 10, height = 7)
  cat("  Saved: metabolomics_KEGG_bubble.pdf (MetaboAnalystR)\n")
}

# KEGGREST-based MSEA (hypergeometric test)
valid_ids <- sig_df_KEGG$KEGG[!is.na(sig_df_KEGG$KEGG)] %>% unique()
cat(sprintf("  Metabolites with KEGG IDs: %d\n", length(valid_ids)))

if (length(valid_ids) >= 3) {
  kegg_pathways <- keggList("pathway", "hsa")
  pathway_ids <- names(kegg_pathways)

  get_pathway_compounds <- function(pid) {
    pw <- tryCatch(keggGet(pid)[[1]], error = function(e) NULL)
    if (is.null(pw) || is.null(pw$COMPOUND)) return(NULL)
    return(names(pw$COMPOUND))
  }

  pathway_compound_list <- lapply(pathway_ids, get_pathway_compounds)
  names(pathway_compound_list) <- pathway_ids
  pathway_compound_list <- pathway_compound_list[!sapply(pathway_compound_list, is.null)]

  all_compounds <- unique(unlist(pathway_compound_list))
  DE_kegg_id <- valid_ids

  msea_list <- lapply(names(pathway_compound_list), function(pid) {
    pw_cmpds <- pathway_compound_list[[pid]]
    overlap <- intersect(DE_kegg_id, pw_cmpds)
    k <- length(overlap)
    M <- length(all_compounds); n <- length(pw_cmpds); N <- length(DE_kegg_id)
    pval <- phyper(q = k - 1, m = n, n = M - n, k = N, lower.tail = FALSE)
    data.frame(Pathway_ID = pid, Pathway_Name = kegg_pathways[pid],
               DE_in_Pathway = k, Pathway_Size = n, Pvalue = pval,
               Overlap_IDs = paste(overlap, collapse = ";"))
  })

  msea_df <- bind_rows(msea_list) %>%
    mutate(FDR = p.adjust(Pvalue, method = "BH"),
           RichFactor = DE_in_Pathway / Pathway_Size) %>% arrange(Pvalue)
  write.csv(msea_df, file.path(OUT_DIR, "MSEA_results.csv"), row.names = FALSE)

  top_msea <- msea_df %>% filter(Pvalue < 0.05) %>%
    slice_max(order_by = RichFactor, n = 20) %>%
    mutate(Pathway_Name = str_replace(Pathway_Name, " - Homo sapiens \\(human\\)", ""),
           Pathway_Name = factor(Pathway_Name, levels = rev(Pathway_Name)))

  if (nrow(top_msea) > 0) {
    p_bubble <- ggplot(top_msea, aes(x = RichFactor, y = Pathway_Name)) +
      geom_point(aes(size = DE_in_Pathway, color = Pvalue)) +
      scale_color_gradient(low = "red", high = "blue", name = "p-value") +
      scale_size(range = c(3, 8), name = "DE count") +
      labs(title = "Metabolite Set Enrichment Analysis (MSEA)",
           x = "Rich Factor", y = "KEGG Pathway") +
      theme_bw(base_size = 14) + xlim(0, max(top_msea$RichFactor) + 0.05)
    ggsave(file.path(OUT_DIR, "MSEA_bubble.pdf"), p_bubble, width = 10, height = 7)
  }
}

# =============================================================================
# PART 7: Transcript-Metabolite Pearson Correlation
# =============================================================================
cat("\n===== PART 7: Transcript-Metabolite Correlation =====\n")

metabolite_names_7 <- c('Beta-Alanine', '1,3-diaminopropane', 'Uracil',
                        'Spermidine', 'Glutamine', "5'-Cytidylic acid",
                        'Dephosphocoenzyme A (Dephospho-CoA)')
exp_sig_select <- exp_sig[intersect(metabolite_names_7, rownames(exp_sig)), ]
if (nrow(exp_sig_select) > 0) {
  exp_sig_select_log <- log10(exp_sig_select + 1)

  gene_names <- c('KMO', 'IL4I1', 'CYP1B1', 'TDO2', 'KYNU', 'TPH1', 'DDC')
  # RNA-seq expression matrix (exp) should be loaded from bulkRNA-seq analysis
  if (exists("exp") && all(gene_names %in% rownames(exp))) {
    exp_gene_select <- exp[gene_names, ]
    colnames(exp_gene_select) <- gsub("^TIGIT", "TIGIT.", gsub("^IgG", "IgG.", colnames(exp_gene_select)))
    combined_exp <- rbind(exp_sig_select_log, exp_gene_select)
    combined_exp <- t(combined_exp) %>% as.data.frame()

    corr_result <- corr.test(combined_exp, method = "pearson", adjust = "fdr", use = "pairwise.complete.obs")
    cor_matrix <- corr_result$r; p_matrix <- corr_result$p
    rownames(cor_matrix) <- colnames(cor_matrix) <- gsub("['()]", "", colnames(combined_exp))
    rownames(p_matrix) <- colnames(p_matrix) <- gsub("['()]", "", colnames(combined_exp))

    pdf(file.path(OUT_DIR, "metabolite_gene_correlation.pdf"), width = 10, height = 8)
    corrplot(cor_matrix, method = "pie", type = "lower", tl.col = "black", tl.srt = 45,
             tl.cex = 0.6, p.mat = p_matrix, sig.level = c(0.001, 0.01, 0.05),
             insig = "label_sig", pch.cex = 1.2, pch.col = "white",
             col = colorRampPalette(c("blue", "white", "red"))(50), mar = c(0, 0, 0, 0))
    dev.off()
    cat("  Saved: metabolite_gene_correlation.pdf\n")
  }
}

# =============================================================================
# PART 8: Pathview — KEGG Pathway Maps (Fig 6D supplement)
# =============================================================================
cat("\n===== PART 8: Pathview KEGG Maps =====\n")

sig_gene_path <- file.path("..", "bulkRNA-seq", "output", "sig_gene.csv")
if (file.exists(sig_gene_path)) {
  sig_gene <- read.csv(sig_gene_path)
  gene_names <- c("KMO", "IL4I1", "CYP1B1", "TDO2", "KYNU", "TPH1", "DDC")
  gene_data <- sig_gene[sig_gene$gene %in% gene_names, ]
  if (nrow(gene_data) > 0) {
    gene_data_matrix <- as.matrix(gene_data$logFC)
    rownames(gene_data_matrix) <- gene_data$ENTREZID

    cpd_data <- sig_df_KEGG[!is.na(sig_df_KEGG$KEGG), ]
    if (nrow(cpd_data) > 0) {
      cpd_data_matrix <- as.matrix(cpd_data$log2FC)
      rownames(cpd_data_matrix) <- cpd_data$KEGG

      pathview(gene.data = gene_data_matrix, cpd.data = cpd_data_matrix,
               pathway.id = "hsa00380", species = "hsa",
               out.suffix = "tryptophan", gene.idtype = "entrez", cpd.idtype = "kegg",
               keys.align = "y", kegg.native = TRUE, same.layer = FALSE,
               low = list(gene = "purple", cpd = "purple"),
               mid = list(gene = "gray", cpd = "gray"),
               high = list(gene = "red", cpd = "red"))
      cat("  Pathview output saved (hsa00380.tryptophan.*)\n")
    }
  }
}

# =============================================================================
# PART 9: Random Forest Biomarker Selection
# =============================================================================
cat("\n===== PART 9: Random Forest =====\n")

meta_exp_rf <- t(meta_exp)
rf_model <- randomForest(meta_exp_rf, y = group_list, importance = TRUE, ntree = 500)
imp_df <- as.data.frame(importance(rf_model, type = 1))
imp_df$Metabolite <- rownames(imp_df)
imp_df <- imp_df[order(imp_df$MeanDecreaseAccuracy, decreasing = TRUE), ]

p_rf <- ggplot(head(imp_df, 20), aes(x = reorder(Metabolite, MeanDecreaseAccuracy),
                                      y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
  labs(title = "Top 20 Important Metabolites (Random Forest)",
       x = "Metabolite", y = "Importance") + theme_minimal(base_size = 14)
ggsave(file.path(OUT_DIR, "RF_top20_metabolites.pdf"), p_rf, width = 10, height = 8)

# =============================================================================
# PART 10: DIABLO — Multi-omics Integration (Fig 6E)
#   mixOmics block.splsda: integrates RNA-seq + metabolomics for joint
#   dimension reduction and feature selection
# =============================================================================
cat("\n===== PART 10: DIABLO Multi-omics Integration (Fig 6E) =====\n")

# Load RNA-seq VST expression from bulkRNA-seq output
exp_vst_path <- file.path("..", "bulkRNA-seq", "output", "exp_vst.csv")
if (file.exists(exp_vst_path)) {
  exp_vst <- as.matrix(read.csv(exp_vst_path, row.names = 1))

  # Align sample columns
  common_samps <- intersect(colnames(exp_vst), colnames(meta_exp))
  if (length(common_samps) >= 4) {
    meta_exp_aligned <- meta_exp[, common_samps, drop = FALSE]
    exp_vst_aligned <- exp_vst[, common_samps, drop = FALSE]

    # Build group factor from column names
    Y <- factor(ifelse(grepl("^IgG", common_samps), "IgG", "TIGIT"),
                levels = c("IgG", "TIGIT"))

    X_rna <- t(exp_vst_aligned)
    X_metabolite <- t(meta_exp_aligned)

    # Top 5000 variable genes
    gene_var <- apply(X_rna, 2, var)
    top_genes <- order(gene_var, decreasing = TRUE)[1:min(5000, ncol(X_rna))]
    X_rna_filtered <- X_rna[, top_genes]

    cat(sprintf("  RNA: %d x %d, Metabolite: %d x %d\n",
        nrow(X_rna_filtered), ncol(X_rna_filtered),
        nrow(X_metabolite), ncol(X_metabolite)))

    # Design matrix (0.1 correlation between blocks)
    design <- matrix(0.1, ncol = 2, nrow = 2,
                     dimnames = list(c("RNA", "Metabolite"), c("RNA", "Metabolite")))
    diag(design) <- 0

    # Tuning (LOO CV)
    cat("  Tuning DIABLO (LOO)...\n")
    tune_res <- tune.block.splsda(
      X = list(RNA = X_rna_filtered, Metabolite = X_metabolite),
      Y = Y, ncomp = 2,
      test.keepX = list(RNA = c(10, 20, 50), Metabolite = c(5, 10, 20)),
      design = design, validation = "loo",
      dist = "centroids.dist", measure = "BER"
    )
    pdf(file.path(OUT_DIR, "DIABLO_tuning.pdf"), width = 8, height = 6)
    plot(tune_res); dev.off()

    # Final model
    final_diablo <- block.splsda(
      X = list(RNA = X_rna_filtered, Metabolite = X_metabolite),
      Y = Y, ncomp = 2, keepX = tune_res$choice.keepX, design = design
    )

    # Export selected features
    select_rna <- selectVar(final_diablo, comp = 1, block = 'RNA')
    write.csv(data.frame(GeneID = select_rna$name, Value = select_rna$value.var),
              file.path(OUT_DIR, "DIABLO_Genes_Top.csv"), row.names = FALSE)
    select_meta <- selectVar(final_diablo, comp = 1, block = 'Metabolite')
    write.csv(data.frame(Metabolite = select_meta$name, Value = select_meta$value.var),
              file.path(OUT_DIR, "DIABLO_Metabolites_Top.csv"), row.names = FALSE)

    # Fig 6E plots
    pdf(file.path(OUT_DIR, "DIABLO_SamplePlot.pdf"), width = 10, height = 8)
    plotIndiv(final_diablo, ind.names = TRUE, ellipse = TRUE, legend = TRUE,
              title = "DIABLO: IgG vs TIGIT", cex = 3)
    dev.off()

    pdf(file.path(OUT_DIR, "DIABLO_LoadingPlot.pdf"), width = 14, height = 8)
    plotLoadings(final_diablo, comp = 1, contrib = 'max', method = 'median',
                 size.name = 1.2, title = "Top Contributors (IgG vs TIGIT)")
    dev.off()

    pdf(file.path(OUT_DIR, "DIABLO_CircosPlot.pdf"), width = 12, height = 12)
    circosPlot(final_diablo, cutoff = 0.9, line = TRUE,
               color.blocks = c('darkorchid', 'brown1'),
               color.cor = c("chocolate3", "grey20"),
               size.labels = 2.0, size.variables = 1.5, size.legend = 1.5)
    dev.off()

    pdf(file.path(OUT_DIR, "DIABLO_Heatmap.pdf"), width = 14, height = 12)
    cimDiablo(final_diablo, color.blocks = c('darkorchid', 'brown1'),
              comp = 1, margin = c(15, 25), legend.position = "right",
              trim = 3, cluster = "none", row.cex = 1.5, col.cex = 1.2)
    dev.off()

    cat("  DIABLO output: SamplePlot, LoadingPlot, CircosPlot, Heatmap\n")
  } else {
    cat(sprintf("  Only %d common samples (need >=4). Skipping DIABLO.\n", length(common_samps)))
  }
} else {
  cat("  RNA-seq VST file not found. Skipping DIABLO.\n")
  cat("  Expected: ", normalizePath(exp_vst_path, mustWork = FALSE), "\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n========================================\n")
cat("  TIGIT Metabolomics + Multi-omics Complete\n")
cat("========================================\n")
cat("  Fig 6A: metabolomics_volcano.pdf\n")
cat("  Fig 6B: MSEA_bubble.pdf / metabolomics_KEGG_bubble.pdf\n")
cat("  Fig 6D: metabolomics_7key_boxplots.pdf + Pathview maps\n")
cat("  Fig 6E: DIABLO_SamplePlot/LoadingPlot/CircosPlot/Heatmap.pdf\n")
cat("  Additional: PCA, heatmaps, RF, metabolite-gene correlation\n")
cat("========================================\n")
