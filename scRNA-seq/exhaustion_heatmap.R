# =============================================================================
# exhaustion_heatmap.R — NK cell exhaustion-associated metabolic pathway
#                        heatmaps (TIGIT High vs Low, 3h + 24h)
# TB-TIGIT-NK manuscript
#
# Parameterized script: set SELECTED_SET below to choose pathway set.
# Available sets:
#   "exhaustion"       — 17 exhaustion-relevant metabolic pathways
#   "mito_energy"      — Mitochondria-related energy metabolism (6 pathways)
#   "mito_energy_redox" — Mito energy + redox (5 pathways)
#   "metabolomics"     — Metabolomics-validated pathways (8 pathways)
# =============================================================================

# ===== SETUP =====
rm(list = ls())
cat("=== Exhaustion Heatmap ===\n")

library(Seurat)
library(readxl)
library(dplyr)
library(pheatmap)

source("utils/scRNA-seq.R")

OUT_DIR <- "./output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
DATA_DIR <- "./data"

set.seed(42)

# ===== SELECT PATHWAY SET =====
# Edit this line to choose which set to analyze:
SELECTED_SET <- "exhaustion"  # options: "exhaustion", "mito_energy", "mito_energy_redox", "metabolomics"

# ===== LOAD DATA =====
cat("\n===== Loading data =====\n")
obj <- readRDS(file.path(DATA_DIR, "mtb_final.rds"))
nk <- obj[, obj$cell_type_lowerres == "NK"]
cat(sprintf("NK cells: %d\n", ncol(nk)))

df <- read_excel(file.path(DATA_DIR, "114\u6761\u4ee3\u8c22\u901a\u8def\u57fa\u56e0.xlsx"),
                 sheet = "Sheet1")

# ===== DEFINE PATHWAY SETS =====
cat("\n===== Defining pathway set:", SELECTED_SET, "=====\n")

pathway_sets <- list(

  # ---------------------------------------------------------------------------
  # Set 1: exhaustion — 17 exhaustion-relevant metabolic pathways
  # ---------------------------------------------------------------------------
  exhaustion = list(
    pathways = list(
      Glycolysis                  = "Glycolysis",
      `Citric Acid Cycle`         = "Citric Acid Cycle",
      `Oxidative Phosphorylation` = "Oxidative Phosphorylation",
      `Pyruvate Metabolism`       = "Pyruvate Metabolism",
      `Fatty Acid Degradation`    = "Fatty Acid Degradation",
      `Tryptophan Metabolism`     = "Tryptophan Metabolism",
      `Kynurenine Metabolism`     = "Kynurenine Metabolism",
      `Arginine Biosynthesis`     = "Arginine Biosynthesis",
      `Arginine-Proline`          = "Arginine and Proline Metabolism",
      `Glutamine-Glutamate`       = "D-Glutamine and D-Glutamate Metabolism",
      `Glutathione Metabolism`    = "Glutathione Metabolism",
      `Pentose Phosphate`         = "Pentose Phosphate",
      `Folate One-Carbon`         = "Folate One Carbon Metabolism",
      `Polyamine Biosynthesis`    = "Polyamine Biosynthesis",
      `Pyrimidine Metabolism`     = "Pyrimidine Metabolism",
      `Purine Metabolism`         = "Purine Metabolism",
      `Methionine Cycle`          = "Methionine Cycle"
    ),
    categories = c(
      Glycolysis                  = "Energy Metabolism",
      `Citric Acid Cycle`         = "Energy Metabolism",
      `Oxidative Phosphorylation` = "Energy Metabolism",
      `Pyruvate Metabolism`       = "Energy Metabolism",
      `Fatty Acid Degradation`    = "Lipid Metabolism",
      `Tryptophan Metabolism`     = "Amino Acid",
      `Kynurenine Metabolism`     = "Amino Acid",
      `Arginine Biosynthesis`     = "Amino Acid",
      `Arginine-Proline`          = "Amino Acid",
      `Glutamine-Glutamate`       = "Amino Acid",
      `Methionine Cycle`          = "Amino Acid",
      `Glutathione Metabolism`    = "Redox",
      `Pentose Phosphate`         = "Nucleotide",
      `Folate One-Carbon`         = "Nucleotide",
      `Polyamine Biosynthesis`    = "Amino Acid",
      `Pyrimidine Metabolism`     = "Nucleotide",
      `Purine Metabolism`         = "Nucleotide"
    ),
    cat_order = c("Energy Metabolism", "Amino Acid", "Lipid Metabolism",
                  "Nucleotide", "Redox"),
    cat_colors = c(
      "Energy Metabolism" = "#E64B35",
      "Amino Acid"        = "#4DBBD5",
      "Lipid Metabolism"  = "#F39B7F",
      "Nucleotide"        = "#00A087",
      "Redox"             = "#3C5488"
    ),
    file_prefix = "NK_TIGIT_Exhaustion_Metabolic",
    plot_title   = "NK Cell Exhaustion-Associated Metabolic Pathways",
    width        = 7,
    height       = NULL,   # auto-calculated
    fontsize_row = 12
  ),

  # ---------------------------------------------------------------------------
  # Set 2: mito_energy — Mitochondria-related energy metabolism (6 pathways)
  # ---------------------------------------------------------------------------
  mito_energy = list(
    pathways = list(
      Glycolysis                  = "Glycolysis",
      `Citric Acid Cycle`         = "Citric Acid Cycle",
      `Oxidative Phosphorylation` = "Oxidative Phosphorylation",
      `Pyruvate Metabolism`       = "Pyruvate Metabolism",
      `Fatty Acid Degradation`    = "Fatty Acid Degradation",
      `Glutathione Metabolism`    = "Glutathione Metabolism"
    ),
    categories = c(
      Glycolysis                  = "Energy Metabolism",
      `Citric Acid Cycle`         = "Energy Metabolism",
      `Oxidative Phosphorylation` = "Energy Metabolism",
      `Pyruvate Metabolism`       = "Energy Metabolism",
      `Fatty Acid Degradation`    = "Lipid Metabolism",
      `Glutathione Metabolism`    = "Redox"
    ),
    cat_order = c("Energy Metabolism", "Lipid Metabolism", "Redox"),
    cat_colors = c(
      "Energy Metabolism" = "#E64B35",
      "Lipid Metabolism"  = "#F39B7F",
      "Redox"             = "#3C5488"
    ),
    file_prefix = "NK_TIGIT_Mito_Energy",
    plot_title   = "NK Cell Mito-Energy Metabolism (TIGIT High vs Low)",
    width        = 7,
    height       = 5,
    fontsize_row = 12
  ),

  # ---------------------------------------------------------------------------
  # Set 3: mito_energy_redox — Mito energy + redox (5 pathways)
  # ---------------------------------------------------------------------------
  mito_energy_redox = list(
    pathways = list(
      Glycolysis                  = "Glycolysis",
      `Pyruvate Metabolism`       = "Pyruvate Metabolism",
      `Citric Acid Cycle`         = "Citric Acid Cycle",
      `Oxidative Phosphorylation` = "Oxidative Phosphorylation",
      `Glutathione Metabolism`    = "Glutathione Metabolism"
    ),
    categories = c(
      Glycolysis                  = "Mitochondrial Energy",
      `Pyruvate Metabolism`       = "Mitochondrial Energy",
      `Citric Acid Cycle`         = "Mitochondrial Energy",
      `Oxidative Phosphorylation` = "Mitochondrial Energy",
      `Glutathione Metabolism`    = "Redox Homeostasis"
    ),
    cat_order = c("Mitochondrial Energy", "Redox Homeostasis"),
    cat_colors = c(
      "Mitochondrial Energy" = "#E64B35",
      "Redox Homeostasis"   = "#3C5488"
    ),
    file_prefix = "NK_TIGIT_Mito_Energy_Redox",
    plot_title   = "NK Cell Mitochondrial Energy & Redox Homeostasis",
    width        = 6.5,
    height       = 4.5,
    fontsize_row = 12
  ),

  # ---------------------------------------------------------------------------
  # Set 4: metabolomics — Metabolomics-validated pathways (8 pathways)
  # ---------------------------------------------------------------------------
  metabolomics = list(
    pathways = list(
      `Glutamine-Glutamate`         = "D-Glutamine and D-Glutamate Metabolism",
      `Alanine-Aspartate-Glutamate` = "Alanine, Aspartate and Glutamate Metabolism",
      `Beta-Alanine Metabolism`     = "Beta-Alanine Metabolism",
      `Arginine-Proline`            = "Arginine and Proline Metabolism",
      `Polyamine Biosynthesis`      = "Polyamine Biosynthesis",
      `Pyrimidine Metabolism`       = "Pyrimidine Metabolism",
      `Pantothenate-CoA`            = "Pantothenate and CoA Biosynthesis",
      `Glutathione Metabolism`      = "Glutathione Metabolism"
    ),
    categories = c(
      `Glutamine-Glutamate`         = "Amino Acid & Derivatives",
      `Alanine-Aspartate-Glutamate` = "Amino Acid & Derivatives",
      `Beta-Alanine Metabolism`     = "Amino Acid & Derivatives",
      `Arginine-Proline`            = "Amino Acid & Derivatives",
      `Polyamine Biosynthesis`      = "Amino Acid & Derivatives",
      `Pyrimidine Metabolism`       = "Nucleotide Metabolism",
      `Pantothenate-CoA`            = "Cofactor & Vitamin",
      `Glutathione Metabolism`      = "Redox Metabolism"
    ),
    cat_order = c("Amino Acid & Derivatives", "Nucleotide Metabolism",
                  "Cofactor & Vitamin", "Redox Metabolism"),
    cat_colors = c(
      "Amino Acid & Derivatives" = "#4DBBD5",
      "Nucleotide Metabolism"    = "#00A087",
      "Cofactor & Vitamin"       = "#E18727",
      "Redox Metabolism"         = "#3C5488"
    ),
    file_prefix = "NK_TIGIT_Metabolomics",
    plot_title   = "NK Cell Metabolomics-Validated Pathways (TIGIT High vs Low)",
    width        = 7.5,
    height       = 5,
    fontsize_row = 11
  )
)

# ===== VALIDATE SELECTION =====
stopifnot(
  "SELECTED_SET must be one of: exhaustion, mito_energy, mito_energy_redox, metabolomics" =
    SELECTED_SET %in% names(pathway_sets)
)

# ===== EXTRACT ACTIVE SET =====
active   <- pathway_sets[[SELECTED_SET]]
pw_list  <- active$pathways
pw_cat   <- active$categories
cat_ord  <- active$cat_order
cat_col  <- active$cat_colors
f_prefix <- active$file_prefix
p_title  <- active$plot_title
p_width  <- active$width
p_height <- active$height
p_fs_row <- active$fontsize_row

# ===== PARSE PATHWAY GENES FROM EXCEL =====
cat("\n===== Parsing pathway genes =====\n")
pathway_genes <- list()
for (nm in names(pw_list)) {
  cn <- pw_list[[nm]]
  if (cn %in% colnames(df)) {
    g <- trimws(na.omit(df[[cn]]))
    g <- g[g != ""]
    g <- intersect(g, rownames(nk))
    if (length(g) >= 2) pathway_genes[[nm]] <- g
  }
}
cat(sprintf("Pathways: %d / %d\n", length(pathway_genes), length(pw_list)))
for (nm in names(pathway_genes)) {
  cat(sprintf("  %s: %d genes\n", nm, length(pathway_genes[[nm]])))
}

# ===== TIGIT GROUPING =====
cat("\n===== TIGIT Low / High group (3h + 24h) =====\n")
r3  <- process_tigit(nk, "3hMTB")
r24 <- process_tigit(nk, "24hMTB")
cat(sprintf("3h:  Low=%d  High=%d\n",
            sum(r3$TIGIT_group == "Low"), sum(r3$TIGIT_group == "High")))
cat(sprintf("24h: Low=%d  High=%d\n",
            sum(r24$TIGIT_group == "Low"), sum(r24$TIGIT_group == "High")))

# ===== GLOBAL Z-SCORES =====
cat("\n===== Computing global z-scores =====\n")
nk_3h_full  <- nk[, nk$timepoint == "3hMTB"]
nk_24h_full <- nk[, nk$timepoint == "24hMTB"]
cat(sprintf("3h full NK:  %d cells\n", ncol(nk_3h_full)))
cat(sprintf("24h full NK: %d cells\n", ncol(nk_24h_full)))

h3  <- compute_pw_global_z(nk_3h_full, r3, pathway_genes)
h24 <- compute_pw_global_z(nk_24h_full, r24, pathway_genes)

common_pw <- intersect(rownames(h3$mat), rownames(h24$mat))
cat(sprintf("Common pathways: %d\n", length(common_pw)))

# ===== CLUSTERING WITHIN CATEGORY (by 24h data) =====
pw_cat_final <- pw_cat[common_pw]
h24_sub <- h24$mat[common_pw, , drop = FALSE]

gene_order <- character()
for (ctg in cat_ord) {
  pw_c <- names(pw_cat_final)[pw_cat_final == ctg]
  if (length(pw_c) > 1) {
    d <- dist(h24_sub[pw_c, , drop = FALSE])
    pw_c <- pw_c[hclust(d, method = "ward.D2")$order]
  }
  if (length(pw_c) > 0) gene_order <- c(gene_order, pw_c)
}

# ===== BUILD COMBINED MATRIX =====
heat_c <- cbind(
  h3$mat[gene_order, c("Low", "High")],
  h24$mat[gene_order, c("Low", "High")]
)
colnames(heat_c) <- c("3h_Low", "3h_High", "24h_Low", "24h_High")
max_abs_c <- max(abs(heat_c), na.rm = TRUE)

# ===== GAPS BETWEEN CATEGORIES =====
pw_cat_ordered <- pw_cat_final[gene_order]
cat_counts <- sapply(cat_ord, function(x) sum(pw_cat_ordered == x))
gap_pos <- cumsum(cat_counts[cat_counts > 0])
if (length(gap_pos) > 1) gap_pos <- gap_pos[-length(gap_pos)]

# ===== ANNOTATIONS =====
ann_col <- data.frame(
  Timepoint = c("3h", "3h", "24h", "24h"),
  TIGIT     = c("Low", "High", "Low", "High"),
  row.names = colnames(heat_c)
)
ann_row <- data.frame(
  Category = factor(pw_cat_ordered, levels = cat_ord),
  row.names = gene_order
)

# ===== DETERMINE HEIGHT =====
if (is.null(p_height)) {
  p_height <- max(4, length(gene_order) * 0.35 + 2)
}

# ===== SAVE HEATMAP PDF =====
cat(sprintf("\n===== Saving heatmap: %s.pdf =====\n", f_prefix))
pdf(file.path(OUT_DIR, paste0(f_prefix, ".pdf")),
    width = p_width, height = p_height)

pheatmap(heat_c,
  color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  breaks = seq(-max_abs_c, max_abs_c, length.out = 101),
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = ann_col, annotation_row = ann_row,
  cellwidth = 24,
  cellheight = max(14, min(20, 200 / length(gene_order))),
  gaps_row = as.numeric(gap_pos), gaps_col = 2,
  annotation_colors = list(
    Timepoint = c("3h" = "#D55E00", "24h" = "#CC79A7"),
    TIGIT     = c(Low = "#E8A838", High = "#C43A3A"),
    Category  = cat_col
  ),
  display_numbers = FALSE,
  fontsize_row = p_fs_row,
  fontsize_col = 12,
  border_color = NA,
  main = p_title
)

dev.off()
cat(sprintf("  >> Saved: %s\n", file.path(OUT_DIR, paste0(f_prefix, ".pdf"))))

# ===== STATS TABLE =====
cat(sprintf("\n===== Saving stats: %s_Stats.csv =====\n", f_prefix))
stats_df <- data.frame(
  Pathway = gene_order,
  H3_Low   = round(heat_c[gene_order, "3h_Low"], 3),
  H3_High  = round(heat_c[gene_order, "3h_High"], 3),
  H24_Low  = round(heat_c[gene_order, "24h_Low"], 3),
  H24_High = round(heat_c[gene_order, "24h_High"], 3),
  Diff_3h  = round(heat_c[gene_order, "3h_High"] - heat_c[gene_order, "3h_Low"], 3),
  Diff_24h = round(heat_c[gene_order, "24h_High"] - heat_c[gene_order, "24h_Low"], 3),
  Category = pw_cat_ordered,
  row.names = NULL
)
stats_df$p_adj_3h  <- round(p.adjust(h3$p[gene_order], "BH"), 4)
stats_df$p_adj_24h <- round(p.adjust(h24$p[gene_order], "BH"), 4)
stats_df$sig <- ifelse(stats_df$p_adj_24h < 0.001, "***",
                ifelse(stats_df$p_adj_24h < 0.01, "**",
                ifelse(stats_df$p_adj_24h < 0.05, "*", "ns")))

print(stats_df, row.names = FALSE)
write.csv(stats_df, file.path(OUT_DIR, paste0(f_prefix, "_Stats.csv")),
          row.names = FALSE)
cat(sprintf("  >> Saved: %s\n",
            file.path(OUT_DIR, paste0(f_prefix, "_Stats.csv"))))

# ===== SUMMARY =====
cat("\n===== Summary =====\n")
n_sig_3h  <- sum(stats_df$p_adj_3h  < 0.05, na.rm = TRUE)
n_sig_24h <- sum(stats_df$p_adj_24h < 0.05, na.rm = TRUE)
n_up_24h  <- sum(stats_df$Diff_24h > 0 & stats_df$p_adj_24h < 0.05, na.rm = TRUE)
n_dn_24h  <- sum(stats_df$Diff_24h < 0 & stats_df$p_adj_24h < 0.05, na.rm = TRUE)
cat(sprintf("  Set: %s\n", SELECTED_SET))
cat(sprintf("  Pathways analyzed: %d\n", nrow(stats_df)))
cat(sprintf("  Significant at 3h:  %d / %d\n", n_sig_3h, nrow(stats_df)))
cat(sprintf("  Significant at 24h: %d / %d\n", n_sig_24h, nrow(stats_df)))
cat(sprintf("  Up in TIGIT High at 24h:  %d\n", n_up_24h))
cat(sprintf("  Down in TIGIT High at 24h: %d\n", n_dn_24h))
cat("\n=== DONE ===\n")
