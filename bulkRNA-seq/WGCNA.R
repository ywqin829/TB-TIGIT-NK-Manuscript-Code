# =============================================================================
# WGCNA.R -- Weighted Gene Co-expression Network Analysis
# TB-TIGIT-NK manuscript
#
# Input:  data/exp.csv       (expression matrix, genes x samples)
#         data/sig_gene.csv  (DEGs with at least gene + logFC columns)
#         data/pd_TCGA.csv   (sample metadata with Group, Age, OS, etc.)
#         data/deg.csv       (full DEG table for GSEA fallback)
#
# Output: (all in ./output/)
#   1. Sample_dendrogram.pdf
#   2. Soft_threshold.pdf (only when recomputed)
#   3. Module-trait_heatmap.pdf
#   4. MM_GS_{turquoise,brown,green}.pdf (module membership scatter plots)
#   5. Network_heatmap.pdf (TOM-based network heatmap)
#   6. GO/KEGG enrichment CSVs + barplot/dotplot PDFs per module
#   7. GSEA results (GO + KEGG)
#
# Dependencies:
#   ../utils/go_kegg.R  (run_go_enrichment, run_kegg_enrichment,
#                        run_gsea, save_enrichment)
# =============================================================================

# ===== SETUP =====
rm(list = ls())
cat("=== WGCNA Analysis ===\n")

library(WGCNA)
library(impute)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(tibble)
library(enrichplot)
library(GOplot)
library(RColorBrewer)
library(gplots)

source("utils/go_kegg.R")

OUT_DIR <- "./output"
DATA_DIR <- "./data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# Allow multi-threading in WGCNA
enableWGCNAThreads()

# ===== LOAD DATA =====
cat("\n===== Loading data =====\n")

exp <- read.csv(file.path(DATA_DIR, "exp.csv"), row.names = 1)
sig_gene <- read.csv(file.path(DATA_DIR, "sig_gene.csv"), row.names = 1)
cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(exp), ncol(exp)))
cat(sprintf("  DEGs: %d genes\n", nrow(sig_gene)))

# Clean gene names: remove rows whose names contain "-" or "-."
non_matching_rows <- !grepl("^-|-\\.", rownames(exp))
exp_cleaned <- exp[non_matching_rows, ]
cat(sprintf("  After removing invalid gene names: %d genes\n", nrow(exp_cleaned)))

# ===== PREPARE EXPRESSION DATA =====
cat("\n===== Preparing expression data =====\n")

# Select top 5000 genes by MAD
mad_order <- order(apply(exp_cleaned, 1, mad), decreasing = TRUE)
datExpr <- t(exp_cleaned[mad_order[1:5000], ])

# Filter out bad samples/genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  cat("  Removing bad samples/genes...\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
cat(sprintf("  datExpr: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

# ===== SAMPLE CLUSTERING =====
cat("\n===== Sample clustering =====\n")

sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file.path(OUT_DIR, "Sample_dendrogram.pdf"), width = 12, height = 6)
par(mar = c(0, 5, 2, 0))
plot(sampleTree,
     main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
cat("  Sample_dendrogram.pdf saved.\n")

# ===== TRAIT DATA =====
cat("\n===== Preparing trait data =====\n")

pdata <- read.csv(file.path(DATA_DIR, "pd_TCGA.csv"), row.names = 1)

# Align rownames between pdata and datExpr
common_rownames <- intersect(rownames(pdata), rownames(datExpr))
pdata <- pdata[common_rownames, , drop = FALSE]
datExpr <- datExpr[common_rownames, , drop = FALSE]

# Ensure identical ordering
pdata <- pdata[order(rownames(pdata)), , drop = FALSE]
datExpr <- datExpr[order(rownames(datExpr)), , drop = FALSE]
stopifnot(identical(rownames(pdata), rownames(datExpr)))
cat(sprintf("  %d samples with trait data\n", nrow(pdata)))

# Build trait matrix (Group + clinical variables)
datTraits <- data.frame(
  IgG       = as.numeric(pdata$Group == "IgG"),
  AlphaTIGIT = as.numeric(pdata$Group == "AlphaTIGIT"),
  Age       = as.numeric(pdata$Age),
  OS        = as.numeric(pdata$OS),
  Weight    = as.numeric(pdata$Weight),
  Status    = as.numeric(factor(pdata$status)),
  T         = as.numeric(factor(pdata$T)),
  Stage     = as.numeric(factor(pdata$Stage)),
  row.names = rownames(pdata)
)

# ===== PART 1: Network Construction (expensive -- pre-compute or load) =====
cat("\n===== PART 1: Network Construction =====\n")

if (FALSE) {
  # --- Soft threshold selection ---
  cat("  Selecting soft threshold...\n")
  powers <- c(1:10, seq(from = 12, to = 30, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  cat(sprintf("  Estimated power: %d\n", sft$powerEstimate))

  pdf(file.path(OUT_DIR, "Soft_threshold.pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  cex1 <- 0.9
  # Scale independence
  plot(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.9, col = "red")
  # Mean connectivity
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       type = "n", main = "Mean connectivity")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5],
       labels = powers, cex = cex1, col = "red")
  dev.off()
  cat("  Soft_threshold.pdf saved.\n")

  # --- blockwiseModules ---
  cat("  Running blockwiseModules (this may take a while)...\n")
  net <- blockwiseModules(
    datExpr,
    power             = sft$powerEstimate,
    TOMType           = "unsigned",
    minModuleSize     = 30,
    reassignThreshold = 0,
    mergeCutHeight    = 0.25,
    deepSplit         = 2,
    numericLabels     = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs          = TRUE,
    saveTOMFileBase   = file.path(OUT_DIR, "testTOM"),
    verbose           = 3
  )

  save(net, datExpr, file = file.path(OUT_DIR, "WGCNA.Rdata"))
  cat("  WGCNA.Rdata saved.\n")
} else {
  cat("  Loading pre-computed WGCNA.Rdata...\n")
  load(file.path(OUT_DIR, "WGCNA.Rdata"))
  cat(sprintf("  Loaded: %d modules\n", length(unique(net$colors))))
}

# ===== PART 2: Module-Trait Correlation =====
cat("\n===== PART 2: Module-Trait Correlation =====\n")

moduleColors <- labels2colors(net$colors)
moduleLabels <- net$colors
nGenes   <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Recompute MEs from colors
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)

# Correlation + p-values
moduleTraitCor    <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Heatmap
pdf(file.path(OUT_DIR, "Module-trait_heatmap.pdf"), width = 10, height = 8)
par(mar = c(4.5, 11, 3, 0.1))
labeledHeatmap(
  Matrix         = moduleTraitCor,
  xLabels        = names(datTraits),
  yLabels        = names(MEs),
  ySymbols       = names(MEs),
  colorLabels    = FALSE,
  colors         = blueWhiteRed(50),
  textMatrix     = textMatrix,
  setStdMargins  = FALSE,
  cex.text       = 0.7,
  zlim           = c(-1, 1),
  main           = "Module-trait relationships"
)
dev.off()
cat("  Module-trait_heatmap.pdf saved.\n")

# ===== PART 3: Module Membership Scatter Plots =====
cat("\n===== PART 3: Module Membership Scatter Plots =====\n")

modNames <- substring(names(MEs), 3)

# Gene module membership (correlation of each gene with each ME)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# Gene significance for the AlphaTIGIT trait (column 2 in datTraits)
instrait <- datTraits[, 2, drop = FALSE]
geneTraitSignificance <- as.data.frame(cor(datExpr, instrait, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(instrait), sep = "")
names(GSPvalue) <- paste("p.GS.", names(instrait), sep = "")

# Scatter plot for modules of interest
target_modules <- c("turquoise", "brown", "green")
for (mod in target_modules) {
  if (!mod %in% modNames) {
    cat(sprintf("  Module '%s' not found, skipping.\n", mod))
    next
  }
  column <- which(names(MEs) == paste("ME", mod, sep = ""))
  moduleGenes <- moduleColors == mod

  pdf(file.path(OUT_DIR, sprintf("MM_GS_%s.pdf", mod)), width = 7, height = 7)
  par(mar = c(5, 5, 4, 2))
  verboseScatterplot(
    abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSignificance[moduleGenes, 1]),
    xlab = paste("Module Membership in", mod, "module"),
    ylab = "Gene significance for AlphaTIGIT",
    main = paste("Module membership vs. gene significance\n", mod),
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
    col = mod
  )
  dev.off()
  cat(sprintf("  MM_GS_%s.pdf saved.\n", mod))
}

# ===== PART 4: Module Gene Extraction =====
cat("\n===== PART 4: Module Gene Extraction =====\n")

# Build per-module gene list from net$colors
merged_colors <- labels2colors(net$colors)
g <- merge(
  as.data.frame(table(net$colors)),
  as.data.frame(table(merged_colors))
)
colnames(g) <- c("Label", "Freq", "Color", "Freq2")
genes_module <- list()
for (i in seq_len(nrow(g))) {
  genes_module[[i]] <- names(net$colors)[net$colors == g$Label[i]]
}
names(genes_module) <- g$Color

# Intersect module genes with DEGs
module_deg <- list()
for (mod in target_modules) {
  if (!mod %in% names(genes_module)) {
    cat(sprintf("  Module '%s' not in gene list, skipping.\n", mod))
    next
  }
  mod_df <- data.frame(gene = genes_module[[mod]], stringsAsFactors = FALSE)
  # sig_gene has rownames = gene; add gene column for join
  sig_gene_df <- sig_gene
  sig_gene_df$gene <- rownames(sig_gene_df)
  module_deg[[mod]] <- inner_join(sig_gene_df, mod_df, by = "gene")
  cat(sprintf("  %s: %d DEGs in module\n", mod, nrow(module_deg[[mod]])))
}

# ===== PART 5: GO/KEGG Enrichment per Module =====
cat("\n===== PART 5: GO/KEGG Enrichment =====\n")

for (mod in target_modules) {
  if (is.null(module_deg[[mod]]) || nrow(module_deg[[mod]]) == 0) {
    cat(sprintf("  Skipping %s (no DEG-module intersection)\n", mod))
    next
  }

  gene_vec <- module_deg[[mod]]$gene
  cat(sprintf("  [%s] %d genes for enrichment\n", mod, length(gene_vec)))

  # --- GO enrichment ---
  cat(sprintf("  [%s] Running GO enrichment...\n", mod))
  ego <- run_go_enrichment(gene_vec)
  save_enrichment(ego, file.path(OUT_DIR, sprintf("GO_enrichment_%s.csv", mod)))

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    pdf(file.path(OUT_DIR, sprintf("GO_barplot_%s.pdf", mod)), width = 10, height = 8)
    print(barplot(ego, showCategory = 15) +
      ggtitle(paste("GO Enrichment -", mod, "Module")))
    dev.off()

    pdf(file.path(OUT_DIR, sprintf("GO_dotplot_%s.pdf", mod)), width = 12, height = 10)
    print(dotplot(ego, split = "ONTOLOGY", showCategory = 6) +
      facet_grid(ONTOLOGY ~ ., scales = "free") +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title  = element_text(size = 14, face = "bold", hjust = 0.5),
        strip.text  = element_text(size = 12)
      ) +
      ggtitle(paste("GO Enrichment -", mod, "Module")))
    dev.off()
    cat(sprintf("  [%s] GO plots saved.\n", mod))
  }

  # --- KEGG enrichment ---
  cat(sprintf("  [%s] Running KEGG enrichment...\n", mod))
  kk <- run_kegg_enrichment(gene_vec)
  save_enrichment(kk, file.path(OUT_DIR, sprintf("KEGG_enrichment_%s.csv", mod)))

  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    pdf(file.path(OUT_DIR, sprintf("KEGG_barplot_%s.pdf", mod)), width = 10, height = 6)
    print(barplot(kk, showCategory = 15,
              title = paste("KEGG Enrichment -", mod, "Module")))
    dev.off()

    # Save readable version (gene symbols instead of ENTREZID)
    kk_read <- tryCatch(
      DOSE::setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) NULL
    )
    save_enrichment(kk_read, file.path(OUT_DIR, sprintf("KEGG_readable_%s.csv", mod)))
    cat(sprintf("  [%s] KEGG plots saved.\n", mod))
  }
}

# ===== PART 6: GSEA =====
cat("\n===== PART 6: GSEA =====\n")

# Build geneList (named logFC vector, sorted decreasing)
if ("logFC" %in% colnames(sig_gene)) {
  deg_df <- sig_gene
  deg_df$gene <- rownames(deg_df)
} else {
  cat("  logFC not in sig_gene, loading deg.csv...\n")
  deg_df <- read.csv(file.path(DATA_DIR, "deg.csv"))
}

# Map SYMBOL to ENTREZID
gene_map <- tryCatch(
  bitr(deg_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
  error = function(e) NULL
)

if (!is.null(gene_map) && nrow(gene_map) > 0) {
  gene_map <- gene_map %>% distinct(SYMBOL, .keep_all = TRUE)
  colnames(gene_map)[colnames(gene_map) == "SYMBOL"] <- "gene"

  temp <- inner_join(gene_map, deg_df, by = "gene")
  temp <- dplyr::select(temp, ENTREZID, logFC)
  geneList <- temp$logFC
  names(geneList) <- as.character(temp$ENTREZID)
  geneList <- geneList[!duplicated(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)
  cat(sprintf("  geneList: %d genes\n", length(geneList)))

  gsea_results <- run_gsea(geneList)

  # GO GSEA
  if (!is.null(gsea_results$go) && nrow(as.data.frame(gsea_results$go)) > 0) {
    save_enrichment(as.data.frame(gsea_results$go),
                    file.path(OUT_DIR, "GSEA_GO.csv"))

    pdf(file.path(OUT_DIR, "GSEA_GO_plots.pdf"), width = 10, height = 8)
    for (j in seq_len(min(4, nrow(gsea_results$go@result)))) {
      p <- gseaplot2(gsea_results$go, geneSetID = j,
                     title = gsea_results$go@result$Description[j],
                     pvalue_table = TRUE)
      print(p)
    }
    dev.off()
    cat("  GSEA GO plots saved.\n")
  }

  # KEGG GSEA
  if (!is.null(gsea_results$kegg) && nrow(as.data.frame(gsea_results$kegg)) > 0) {
    save_enrichment(as.data.frame(gsea_results$kegg),
                    file.path(OUT_DIR, "GSEA_KEGG.csv"))

    pdf(file.path(OUT_DIR, "GSEA_KEGG_plots.pdf"), width = 10, height = 8)
    for (j in seq_len(min(4, nrow(gsea_results$kegg@result)))) {
      p <- gseaplot2(gsea_results$kegg, geneSetID = j,
                     title = gsea_results$kegg@result$Description[j],
                     pvalue_table = TRUE)
      print(p)
    }
    dev.off()
    cat("  GSEA KEGG plots saved.\n")
  }
} else {
  cat("  WARNING: gene mapping failed, skipping GSEA.\n")
}

# ===== PART 7: TOM Network Heatmap =====
cat("\n===== PART 7: TOM Network Heatmap =====\n")

nSelect <- 400
set.seed(10)
select <- sample(nGenes, size = min(nSelect, nGenes))

dissTOM <- 1 - TOMsimilarityFromExpr(datExpr, power = 6)
selectTOM <- dissTOM[select, select]
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[select]

myheatcol <- colorpanel(250, "red", "orange", "lemonchiffon")
plotDiss <- selectTOM^7
diag(plotDiss) <- NA

pdf(file.path(OUT_DIR, "Network_heatmap.pdf"), width = 10, height = 10)
TOMplot(plotDiss, selectTree, selectColors, col = myheatcol,
        main = "Network heatmap plot, selected genes")
dev.off()
cat("  Network_heatmap.pdf saved.\n")

# ===== SUMMARY =====
cat("\n===== Files saved in ", OUT_DIR, " =====\n")
cat("  1. Sample_dendrogram.pdf\n")
cat("  2. Soft_threshold.pdf (only when recomputed)\n")
cat("  3. Module-trait_heatmap.pdf\n")
cat("  4. MM_GS_{turquoise,brown,green}.pdf\n")
cat("  5. GO_enrichment_{turquoise,brown,green}.csv\n")
cat("  6. GO_barplot_{turquoise,brown,green}.pdf\n")
cat("  7. GO_dotplot_{turquoise,brown,green}.pdf\n")
cat("  8. KEGG_enrichment_{turquoise,brown,green}.csv\n")
cat("  9. KEGG_barplot_{turquoise,brown,green}.pdf\n")
cat(" 10. KEGG_readable_{turquoise,brown,green}.csv\n")
cat(" 11. GSEA_GO.csv / GSEA_GO_plots.pdf\n")
cat(" 12. GSEA_KEGG.csv / GSEA_KEGG_plots.pdf\n")
cat(" 13. Network_heatmap.pdf\n")
cat("=== WGCNA Complete ===\n")
