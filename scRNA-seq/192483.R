# =============================================================================
# 192483.R -- GSE192483 + GSE227691 scRNA-seq Pre-processing and T-cell DEG
# TB-TIGIT-NK manuscript
#
# Datasets:
#   - GSE192483: Lung tissue from TB patients (High_FDG vs Low_FDG)
#   - GSE227691: Additional validation dataset (optional)
#
# Panels:
#   1. QC + Harmony integration
#   2. Cell-type annotation (DotPlot + tSNE)
#   3. Visualization (gene expression, heatmap)
#   4. T-cell subset extraction and DEG (High_FDG vs Low_FDG)
#   5. Activation/Exhaustion module scores
# =============================================================================

# ===== SETUP =====
rm(list = ls())
cat("=== GSE192483 Analysis ===\n")

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(ggpubr)
library(clustree)
library(harmony)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GEOquery)

# Output directory
OUT_DIR <- "./output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("  Output directory:", OUT_DIR, "\n")

# Data directory (relative path -- user should set symlink or copy data here)
DATA_DIR <- "./data/GSE192483_RAW"
stopifnot("Data directory not found" = dir.exists(DATA_DIR))

# Set seed for reproducibility
set.seed(42)

# =============================================================================
# PART 1: QC + Harmony Integration (one-time compute)
# =============================================================================
cat("\n===== PART 1: QC + Harmony Integration =====\n")

if (FALSE) {
  # --- 1a. Organize files ---
  cat("  [1a] Organizing 10X files...\n")
  file_list <- list.files(DATA_DIR,
    pattern = '^GSM.*(barcodes|features|matrix).*\\.gz$', full.names = FALSE)
  sample_names <- str_extract(file_list, "^GSM\\d+_\\w+") %>% unique()

  lapply(sample_names, function(x) {
    y <- file_list[grepl(x, file_list)]
    folder <- file.path(DATA_DIR, x)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    barcode_file <- y[grepl("barcodes", y)]
    features_file <- y[grepl("features", y)]
    matrix_file <- y[grepl("matrix.mtx", y)]
    file.rename(file.path(DATA_DIR, barcode_file), file.path(folder, "barcodes.tsv.gz"))
    file.rename(file.path(DATA_DIR, features_file), file.path(folder, "features.tsv.gz"))
    file.rename(file.path(DATA_DIR, matrix_file), file.path(folder, "matrix.mtx.gz"))
  })

  # --- 1b. Read samples ---
  cat("  [1b] Reading 10X samples...\n")
  sample_dirs <- list.dirs(path = DATA_DIR, full.names = TRUE, recursive = FALSE)
  sceList <- lapply(sample_dirs, function(folder) {
    sample_id <- basename(folder)
    cat("    Reading:", sample_id, "\n")
    mat <- Read10X(folder)
    CreateSeuratObject(counts = mat, project = sample_id, min.cells = 5, min.features = 300)
  })
  names(sceList) <- basename(sample_dirs)

  # --- 1c. Merge ---
  cat("  [1c] Merging samples...\n")
  sce.all <- merge(x = sceList[[1]], y = sceList[-1], add.cell.ids = names(sceList))
  sce.all <- JoinLayers(sce.all)
  cat("    Total cells:", ncol(sce.all), "\n")

  # --- 1d. Metadata from GEO ---
  cat("  [1d] Fetching GEO metadata...\n")
  gse <- getGEO("GSE192483", GSEMatrix = TRUE)[[1]]
  meta_data <- pData(gse) %>%
    select(geo_accession, title, characteristics_ch1) %>%
    rename(
      Sample_geo_accession = geo_accession,
      Sample_title = title,
      Sample_characteristics_ch1 = characteristics_ch1
    ) %>%
    mutate(
      Group = str_extract(Sample_title, "H|L"),
      Patient = str_extract(Sample_title, "^SP\\d+"),
      Tissue_Type = case_when(
        Sample_characteristics_ch1 == "tissue: lung tissue with 18F-FDG avidity" ~ "High_FDG",
        Sample_characteristics_ch1 == "tissue: lung tissue" ~ "Low_FDG",
        TRUE ~ "Unknown"
      )
    )
  cat("    Metadata groups:", table(meta_data$Tissue_Type), "\n")

  sce.all$Sample_geo_accession <- str_extract(sce.all$orig.ident, "^GSM\\d+")
  merged_meta <- merge(sce.all@meta.data, meta_data,
    by = "Sample_geo_accession", all.x = TRUE)
  rownames(merged_meta) <- rownames(sce.all@meta.data)
  sce.all@meta.data <- merged_meta

  saveRDS(sce.all, file.path(OUT_DIR, "sce_GSE192483_raw.rds"))
  cat("  [1d] Saved: sce_GSE192483_raw.rds\n")

  # --- 1e. QC metrics ---
  cat("  [1e] Calculating QC metrics...\n")
  mito_genes <- grep("^MT-", rownames(sce.all), value = TRUE, ignore.case = TRUE)
  sce.all <- PercentageFeatureSet(sce.all, features = mito_genes, col.name = "percent_mito")
  ribo_genes <- grep("^RP[SL]", rownames(sce.all), value = TRUE, ignore.case = TRUE)
  sce.all <- PercentageFeatureSet(sce.all, features = ribo_genes, col.name = "percent_ribo")
  hb_genes <- grep("^HB[A-Z]", rownames(sce.all), value = TRUE, ignore.case = TRUE)
  sce.all <- PercentageFeatureSet(sce.all, features = hb_genes, col.name = "percent_hb")

  p1 <- VlnPlot(sce.all, group.by = "orig.ident",
    features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb"),
    pt.size = 0, ncol = 3) + NoLegend()
  ggsave(file.path(OUT_DIR, "Vlnplot_QC_prefilter.pdf"), p1, width = 16, height = 8)

  # --- 1f. Filter ---
  cat("  [1f] Filtering cells...\n")
  sce.all.filt <- subset(sce.all,
    subset = nFeature_RNA > 200 & nCount_RNA > 500 &
      percent_mito < 25 & percent_ribo > 3 & percent_hb < 1)
  cat("    Cells before filtering:", ncol(sce.all), "\n")
  cat("    Cells after filtering:", ncol(sce.all.filt), "\n")

  # --- 1g. Normalize + HVG + Scale + PCA + Harmony ---
  cat("  [1g] Normalize, HVG, Scale, PCA, Harmony...\n")
  sce.all.filt <- NormalizeData(sce.all.filt)
  sce.all.filt <- FindVariableFeatures(sce.all.filt, nfeatures = 2000)
  sce.all.filt <- ScaleData(sce.all.filt)
  sce.all.filt <- RunPCA(sce.all.filt)
  sce.all.filt <- RunHarmony(sce.all.filt, group.by.vars = "orig.ident")
  cat("    Harmony integration complete.\n")

  # --- 1h. UMAP + tSNE + Clustering ---
  cat("  [1h] UMAP, tSNE, Clustering...\n")
  sce.all.filt <- RunUMAP(sce.all.filt, reduction = "harmony", dims = 1:10)
  sce.all.filt <- RunTSNE(sce.all.filt, reduction = "harmony", dims = 1:10, perplexity = 30)

  for (res in c(0.1, 0.3, 0.5, 0.8, 1)) {
    sce.all.filt <- FindClusters(sce.all.filt, resolution = res)
  }
  sce.all.filt <- FindClusters(sce.all.filt, resolution = 0.4)

  # --- 1i. Save ---
  saveRDS(sce.all.filt, file.path(OUT_DIR, "sce_GSE192483_filtered_harmony.rds"))
  cat("  [1i] Saved: sce_GSE192483_filtered_harmony.rds\n")

  # --- 1j. Clustering tree ---
  cat("  [1j] Generating clustree...\n")
  pdf(file.path(OUT_DIR, "Clustree.pdf"), width = 10, height = 8)
  print(clustree(sce.all.filt, prefix = "RNA_snn_res."))
  dev.off()
  cat("  [1j] Saved: Clustree.pdf\n")

  cat("  === Part 1 complete ===\n")
} else {
  cat("  Skipping Part 1 (QC+Harmony). Set to TRUE to rerun.\n")
}

# =============================================================================
# PART 2: Cell-Type Annotation (one-time compute)
# =============================================================================
cat("\n===== PART 2: Cell-Type Annotation =====\n")

if (FALSE) {
  cat("  [2a] Loading filtered object...\n")
  sce <- readRDS(file.path(OUT_DIR, "sce_GSE192483_filtered_harmony.rds"))
  cat("    Cells:", ncol(sce), "| Genes:", nrow(sce), "\n")

  # --- 2b. Define marker genes per lineage ---
  cat("  [2b] Defining marker gene sets...\n")
  markers <- list(
    Epithelial = c("EPCAM", "KRT19", "KRT18"),
    Alveolar_Type1 = c("AGER", "PDPN", "CLIC5"),
    Alveolar_Type2 = c("SFTPC", "SFTPB", "SFTPA1", "LPCAT1"),
    Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"),
    Fibroblast = c("COL1A1", "DCN", "LUM", "COL1A2"),
    Smooth_Muscle = c("ACTA2", "MYH11", "TAGLN"),
    Myeloid = c("CD14", "CD68", "AIF1", "FCGR3A"),
    Macrophage = c("CD163", "MRC1", "C1QA", "C1QB"),
    Monocyte = c("FCN1", "S100A8", "S100A9", "LYZ"),
    Dendritic = c("CLEC10A", "FCER1A", "CD1C"),
    pDC = c("LILRA4", "CLEC4C", "IL3RA"),
    T_Cell = c("CD3D", "CD3E", "CD2", "CD7"),
    CD4_T = c("CD3D", "CD4", "IL7R"),
    CD8_T = c("CD3D", "CD8A", "CD8B", "GZMK"),
    NK = c("NKG7", "GNLY", "KLRD1", "KLRF1"),
    B_Cell = c("CD79A", "MS4A1", "CD19"),
    Plasma = c("MZB1", "SDC1", "JCHAIN"),
    Neutrophil = c("FCGR3B", "CSF3R", "S100A12"),
    Mast = c("KIT", "TPSAB1", "CPA3"),
    Proliferating = c("MKI67", "TOP2A", "STMN1")
  )

  # --- 2c. DotPlot of markers across clusters ---
  cat("  [2c] Generating marker DotPlot...\n")
  all_marker_genes <- unique(unlist(markers))
  p_dot <- DotPlot(sce, features = all_marker_genes, assay = "RNA",
    group.by = "RNA_snn_res.0.4") +
    RotatedAxis() +
    labs(title = "Marker Gene Expression by Cluster (res 0.4)") +
    theme(axis.text.x = element_text(size = 8))
  ggsave(file.path(OUT_DIR, "DotPlot_markers.pdf"), p_dot, width = 24, height = 6)

  # --- 2d. Module score-based annotation ---
  cat("  [2d] Computing module scores...\n")
  for (cell_type in names(markers)) {
    genes_present <- intersect(markers[[cell_type]], rownames(sce))
    if (length(genes_present) > 1) {
      sce <- AddModuleScore(sce, features = list(genes_present),
        name = paste0("Score_", cell_type), ctrl = 50)
    } else {
      cat("    Warning: <2 genes for", cell_type, "--skipping\n")
    }
  }

  # --- 2e. Assign cell types based on max module score ---
  cat("  [2e] Assigning cell types...\n")
  score_cols <- grep("^Score_", colnames(sce@meta.data), value = TRUE)
  score_matrix <- as.matrix(sce@meta.data[, score_cols])
  colnames(score_matrix) <- gsub("^Score_", "", colnames(score_matrix))
  max_score_type <- apply(score_matrix, 1, function(x) names(which.max(x)))
  sce$CellType <- max_score_type
  cat("    Cell type distribution:\n")
  print(table(sce$CellType))

  # --- 2f. tSNE colored by cell type ---
  cat("  [2f] Generating tSNE by cell type...\n")
  p_tsne <- DimPlot(sce, reduction = "tsne", group.by = "CellType",
    label = TRUE, pt.size = 0.3, repel = TRUE) +
    labs(title = "tSNE -- Cell Type Annotation") +
    theme(legend.text = element_text(size = 10))
  ggsave(file.path(OUT_DIR, "tSNE_celltype.pdf"), p_tsne, width = 10, height = 8)

  p_tsne_cluster <- DimPlot(sce, reduction = "tsne", group.by = "RNA_snn_res.0.4",
    label = TRUE, pt.size = 0.3, repel = TRUE) +
    labs(title = "tSNE -- Clusters (res 0.4)")
  ggsave(file.path(OUT_DIR, "tSNE_cluster.pdf"), p_tsne_cluster, width = 10, height = 8)

  # --- 2g. Save annotated object ---
  saveRDS(sce, file.path(OUT_DIR, "sce_GSE192483_annotated.rds"))
  cat("  [2g] Saved: sce_GSE192483_annotated.rds\n")

  # --- 2h. Export annotation table ---
  anno_summary <- sce@meta.data %>%
    group_by(RNA_snn_res.0.4, CellType) %>%
    summarise(CellCount = n(), .groups = "drop") %>%
    group_by(RNA_snn_res.0.4) %>%
    mutate(Proportion = CellCount / sum(CellCount)) %>%
    arrange(RNA_snn_res.0.4, desc(Proportion))
  write.csv(anno_summary, file.path(OUT_DIR, "cluster_celltype_summary.csv"),
    row.names = FALSE)
  cat("  [2h] Saved: cluster_celltype_summary.csv\n")

  cat("  === Part 2 complete ===\n")
} else {
  cat("  Skipping Part 2 (Annotation). Set to TRUE to rerun.\n")
}

# =============================================================================
# PART 3: Visualization (active by default)
# =============================================================================
cat("\n===== PART 3: Visualization =====\n")

# Load annotated object (must exist from Part 2)
if (file.exists(file.path(OUT_DIR, "sce_GSE192483_annotated.rds"))) {
  cat("  [3a] Loading annotated object...\n")
  sce <- readRDS(file.path(OUT_DIR, "sce_GSE192483_annotated.rds"))
  cat("    Cells:", ncol(sce), "| Genes:", nrow(sce), "\n")

  # --- 3b. tSNE split by Tissue_Type ---
  cat("  [3b] Generating tSNE by Tissue_Type...\n")
  if ("Tissue_Type" %in% colnames(sce@meta.data)) {
    p_tsne_split <- DimPlot(sce, reduction = "tsne", group.by = "CellType",
      split.by = "Tissue_Type", pt.size = 0.3, ncol = 2) +
      labs(title = "tSNE -- High_FDG vs Low_FDG")
    ggsave(file.path(OUT_DIR, "tSNE_split_tissue.pdf"),
      p_tsne_split, width = 16, height = 7)
  }

  # --- 3c. Key gene expression tSNE ---
  cat("  [3c] Generating gene expression tSNE plots...\n")
  key_genes <- c("CD3D", "CD8A", "CD4", "NKG7", "GNLY",
    "CD68", "CD14", "EPCAM", "PECAM1", "MKI67")
  genes_present <- intersect(key_genes, rownames(sce))
  p_genes <- FeaturePlot(sce, features = genes_present, reduction = "tsne",
    pt.size = 0.2, ncol = 4) &
    scale_color_viridis_c(option = "C") &
    theme(plot.title = element_text(size = 10))
  ggsave(file.path(OUT_DIR, "tSNE_key_genes.pdf"), p_genes, width = 18, height = 14)

  # --- 3d. Cell type proportions by condition ---
  cat("  [3d] Cell type proportions by condition...\n")
  if ("Tissue_Type" %in% colnames(sce@meta.data)) {
    prop_df <- sce@meta.data %>%
      group_by(Tissue_Type, CellType) %>%
      summarise(CellCount = n(), .groups = "drop") %>%
      group_by(Tissue_Type) %>%
      mutate(Proportion = CellCount / sum(CellCount) * 100)

    p_prop <- ggplot(prop_df, aes(x = Tissue_Type, y = Proportion, fill = CellType)) +
      geom_bar(stat = "identity", position = "fill", width = 0.6) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(title = "Cell Type Proportions", y = "Proportion", x = "") +
      theme_classic() +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(color = "black", size = 12, face = "bold"),
        axis.text.y = element_text(color = "black", size = 12, face = "bold"),
        axis.title = element_text(face = "bold", size = 14),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 10)
      )
    ggsave(file.path(OUT_DIR, "Celltype_proportions.pdf"), p_prop, width = 8, height = 6)
  }

  # --- 3e. DoHeatmap of marker expression ---
  cat("  [3e] Generating DoHeatmap...\n")
  top_markers <- c(
    "EPCAM", "KRT19",           # Epithelial
    "AGER", "PDPN",             # AT1
    "SFTPC", "SFTPB",           # AT2
    "PECAM1", "VWF",            # Endothelial
    "COL1A1", "DCN",            # Fibroblast
    "CD68", "CD14",             # Myeloid
    "CD163", "C1QA",            # Macrophage
    "FCN1", "S100A8",           # Monocyte
    "CLEC10A", "FCER1A",        # DC
    "CD3D", "CD2",              # T cell
    "CD4", "IL7R",              # CD4 T
    "CD8A", "CD8B",             # CD8 T
    "NKG7", "GNLY",             # NK
    "CD79A", "MS4A1",           # B
    "MZB1", "JCHAIN",           # Plasma
    "FCGR3B", "CSF3R",          # Neutrophil
    "KIT", "TPSAB1"             # Mast
  )
  top_markers <- intersect(top_markers, rownames(sce))
  if (length(top_markers) > 2) {
    p_heat <- DoHeatmap(sce, features = top_markers, group.by = "CellType",
      size = 3, assay = "RNA", slot = "data") +
      scale_fill_gradientn(colors = c("navy", "white", "firebrick")) +
      labs(title = "Marker Gene Heatmap by Cell Type")
    ggsave(file.path(OUT_DIR, "Heatmap_markers.pdf"), p_heat, width = 10, height = 12)
    cat("  [3e] Saved: Heatmap_markers.pdf\n")
  }

  cat("  === Part 3 complete ===\n")
} else {
  cat("  [3] Skipping: annotated object not found. Run Part 2 first.\n")
}

# =============================================================================
# PART 4: T-cell Subset DEG (active by default)
# =============================================================================
cat("\n===== PART 4: T-cell Subset DEG =====\n")

if (file.exists(file.path(OUT_DIR, "sce_GSE192483_annotated.rds"))) {
  cat("  [4a] Loading annotated object...\n")
  sce <- readRDS(file.path(OUT_DIR, "sce_GSE192483_annotated.rds"))

  # --- 4b. Subset T cells ---
  cat("  [4b] Subsetting T cells...\n")
  t_cell_types <- c("T_Cell", "CD4_T", "CD8_T")
  valid_types <- intersect(t_cell_types, unique(sce$CellType))
  cat("    Subsetting CellType in:", paste(valid_types, collapse = ", "), "\n")
  sce_t <- subset(sce, subset = CellType %in% valid_types)
  cat("    T cells extracted:", ncol(sce_t), "\n")

  # --- 4c. Re-normalize T-cell subset ---
  cat("  [4c] Re-normalizing T-cell subset...\n")
  sce_t <- NormalizeData(sce_t, verbose = FALSE)
  sce_t <- FindVariableFeatures(sce_t, nfeatures = 1500, verbose = FALSE)
  sce_t <- ScaleData(sce_t, verbose = FALSE)

  # --- 4d. Check condition balance ---
  cat("  [4d] Condition balance:\n")
  if ("Tissue_Type" %in% colnames(sce_t@meta.data)) {
    print(table(sce_t$Tissue_Type))
  } else {
    cat("    Warning: Tissue_Type not in metadata. Check orig.ident:\n")
    print(table(sce_t$orig.ident))
  }

  # --- 4e. DEG: High_FDG vs Low_FDG ---
  cat("  [4e] Running FindMarkers (High_FDG vs Low_FDG)...\n")
  if ("Tissue_Type" %in% colnames(sce_t@meta.data)) {
    Idents(sce_t) <- "Tissue_Type"
    deg <- FindMarkers(sce_t,
      ident.1 = "High_FDG",
      ident.2 = "Low_FDG",
      min.pct = 0.1,
      logfc.threshold = 0.1,
      verbose = TRUE)
    cat("    DEGs found:", nrow(deg), "\n")
    deg$gene <- rownames(deg)
    deg <- deg %>% arrange(p_val_adj)
    write.csv(deg, file.path(OUT_DIR, "DEG_Tcell_HighFDG_vs_LowFDG.csv"),
      row.names = FALSE)
    cat("  [4e] Saved: DEG_Tcell_HighFDG_vs_LowFDG.csv\n")

    # --- 4f. Volcano plot ---
    cat("  [4f] Generating Volcano plot...\n")
    deg_for_volcano <- deg %>%
      mutate(
        Significance = case_when(
          p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Up",
          p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down",
          TRUE ~ "NS"
        )
      )
    cat("    Up:", sum(deg_for_volcano$Significance == "Up"), "\n")
    cat("    Down:", sum(deg_for_volcano$Significance == "Down"), "\n")

    # Select top genes to label
    top_up <- deg_for_volcano %>%
      filter(Significance == "Up") %>%
      slice_min(p_val_adj, n = 15) %>%
      pull(gene)
    top_down <- deg_for_volcano %>%
      filter(Significance == "Down") %>%
      slice_min(p_val_adj, n = 15) %>%
      pull(gene)
    label_genes <- c(top_up, top_down)

    p_volcano <- EnhancedVolcano(deg_for_volcano,
      lab = deg_for_volcano$gene,
      selectLab = label_genes,
      x = "avg_log2FC",
      y = "p_val_adj",
      pCutoff = 0.05,
      FCcutoff = 0.25,
      pointSize = 2,
      labSize = 4,
      col = c("grey70", "royalblue", "firebrick", "purple"),
      colAlpha = 0.7,
      legendPosition = "top",
      legendLabSize = 10,
      legendIconSize = 4,
      drawConnectors = TRUE,
      widthConnectors = 0.3,
      title = "T-cell DEG: High_FDG vs Low_FDG",
      subtitle = paste0(ncol(sce_t), " T cells | p_adj < 0.05, |log2FC| > 0.25"),
      caption = paste0("Total DEGs: ", nrow(deg)))
    ggsave(file.path(OUT_DIR, "Volcano_Tcell_DEG.pdf"), p_volcano,
      width = 10, height = 10)
    cat("  [4f] Saved: Volcano_Tcell_DEG.pdf\n")

    # --- 4g. Top DEG expression violin ---
    cat("  [4g] Top DEG expression violin plots...\n")
    top_genes <- deg %>%
      filter(p_val_adj < 0.05) %>%
      slice_max(abs(avg_log2FC), n = 12) %>%
      pull(gene)
    if (length(top_genes) > 0) {
      top_genes <- intersect(top_genes, rownames(sce_t))
      if (length(top_genes) > 0) {
        p_violin <- VlnPlot(sce_t, features = top_genes, group.by = "Tissue_Type",
          pt.size = 0, ncol = 4) & NoLegend()
        ggsave(file.path(OUT_DIR, "Vln_TopDEGs.pdf"), p_violin, width = 16, height = 10)
        cat("  [4g] Saved: Vln_TopDEGs.pdf\n")
      }
    }

    # --- 4h. Save T-cell subset ---
    saveRDS(sce_t, file.path(OUT_DIR, "sce_GSE192483_Tcells.rds"))
    cat("  [4h] Saved: sce_GSE192483_Tcells.rds\n")
  } else {
    cat("  [4e] Skipping DEG: Tissue_Type not found in metadata.\n")
  }

  cat("  === Part 4 complete ===\n")
} else {
  cat("  [4] Skipping: annotated object not found. Run Part 2 first.\n")
}

# =============================================================================
# PART 5: Activation/Exhaustion Module Scores (active by default)
# =============================================================================
cat("\n===== PART 5: Activation/Exhaustion Module Scores =====\n")

if (file.exists(file.path(OUT_DIR, "sce_GSE192483_annotated.rds"))) {
  cat("  [5a] Loading annotated object...\n")
  sce <- readRDS(file.path(OUT_DIR, "sce_GSE192483_annotated.rds"))

  # --- 5b. Define gene modules ---
  cat("  [5b] Defining gene modules...\n")
  activation_genes <- c("CD69", "CD25", "IL2RA", "HLA-DRA", "HLA-DRB1",
    "CD38", "CD28", "ICOS", "TNFRSF4", "TNFRSF9")
  exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2",
    "TOX", "TOX2", "ENTPD1", "CD244", "EOMES")
  cytotoxicity_genes <- c("GZMA", "GZMB", "GZMK", "GZMH", "PRF1",
    "GNLY", "NKG7", "FASLG", "TNF", "IFNG")
  proliferation_genes <- c("MKI67", "TOP2A", "STMN1", "PCNA",
    "TYMS", "RRM2", "CDK1", "BIRC5")

  modules <- list(
    Activation = activation_genes,
    Exhaustion = exhaustion_genes,
    Cytotoxicity = cytotoxicity_genes,
    Proliferation = proliferation_genes
  )

  # --- 5c. Compute module scores ---
  cat("  [5c] Computing module scores...\n")
  for (mod_name in names(modules)) {
    genes_present <- intersect(modules[[mod_name]], rownames(sce))
    if (length(genes_present) > 1) {
      sce <- AddModuleScore(sce, features = list(genes_present),
        name = paste0("Module_", mod_name, "_"), ctrl = 100)
      cat("    Added module:", mod_name, "(", length(genes_present), "genes)\n")
    } else {
      cat("    Warning: <2 genes for module", mod_name, "--skipping\n")
    }
  }

  # --- 5d. Violin plots by CellType and Condition ---
  cat("  [5d] Module score violin plots...\n")
  module_cols <- grep("^Module_", colnames(sce@meta.data), value = TRUE)
  # Clean module names (AddModuleScore appends "1" to name)
  names(module_cols) <- gsub("^Module_|_1$", "", module_cols)

  for (mod_name in names(module_cols)) {
    col <- module_cols[mod_name]
    # Rename for clean plotting
    sce[[paste0("Score_", mod_name)]] <- sce[[col]]
  }
  score_cols_clean <- grep("^Score_", colnames(sce@meta.data), value = TRUE)

  # Violin: by CellType
  for (score_name in score_cols_clean) {
    display_name <- gsub("^Score_", "", score_name)
    p_violin <- VlnPlot(sce, features = score_name, group.by = "CellType",
      pt.size = 0) & NoLegend() &
      labs(title = paste0(display_name, " Score by Cell Type"), y = "Score")
    ggsave(file.path(OUT_DIR, paste0("Vln_", display_name, "_byCellType.pdf")),
      p_violin, width = 10, height = 5)
  }
  cat("    Saved: Violin plots by CellType\n")

  # Violin: by Tissue_Type (T cells only)
  if ("Tissue_Type" %in% colnames(sce@meta.data)) {
    t_cell_types <- c("T_Cell", "CD4_T", "CD8_T")
    valid_types <- intersect(t_cell_types, unique(sce$CellType))
    if (length(valid_types) > 0) {
      sce_t <- subset(sce, subset = CellType %in% valid_types)
      cat("    T-cell subset for module comparison:", ncol(sce_t), "cells\n")

      for (score_name in score_cols_clean) {
        display_name <- gsub("^Score_", "", score_name)
        # Add score as a named column for plotting
        sce_t[[display_name]] <- sce_t[[score_name]]
        p_box <- ggboxplot(sce_t@meta.data,
          x = "Tissue_Type", y = display_name,
          fill = "Tissue_Type", palette = c("#E64B35", "#4DBBD5"),
          add = "jitter", add.size = 0.3, size = 0.5) +
          stat_compare_means(method = "wilcox.test", label = "p.signif",
            hide.ns = TRUE, size = 5) +
          labs(title = paste0(display_name, " Score (T cells)"),
            y = "Module Score", x = "") +
          theme_classic() +
          theme(
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"),
            axis.text.x = element_text(color = "black", size = 12, face = "bold"),
            axis.text.y = element_text(color = "black", size = 12, face = "bold"),
            axis.title = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
            legend.position = "none"
          )
        ggsave(file.path(OUT_DIR,
          paste0("Boxplot_", display_name, "_Tcell_HighFDGvsLowFDG.pdf")),
          p_box, width = 5, height = 6)
      }
      cat("    Saved: Boxplots by Tissue_Type (T cells)\n")
    }
  }

  # --- 5e. Feature plots on tSNE ---
  cat("  [5e] Module score feature plots (tSNE)...\n")
  score_features <- score_cols_clean[score_cols_clean %in% colnames(sce@meta.data)]
  if (length(score_features) > 0 && length(score_features) <= 9) {
    p_features <- FeaturePlot(sce, features = score_features, reduction = "tsne",
      pt.size = 0.2, ncol = 3) &
      scale_color_viridis_c(option = "D") &
      theme(plot.title = element_text(size = 11))
    ggsave(file.path(OUT_DIR, "tSNE_module_scores.pdf"), p_features,
      width = 15, height = 4 * ceiling(length(score_features) / 3))
    cat("    Saved: tSNE_module_scores.pdf\n")
  } else if (length(score_features) > 9) {
    # Save separately
    n_plots <- length(score_features)
    n_col <- 3
    n_row <- ceiling(n_plots / n_col)
    for (page in seq_len(ceiling(n_row / 3))) {
      idx_start <- (page - 1) * 9 + 1
      idx_end <- min(page * 9, n_plots)
      feat_subset <- score_features[idx_start:idx_end]
      p_features <- FeaturePlot(sce, features = feat_subset, reduction = "tsne",
        pt.size = 0.2, ncol = 3) &
        scale_color_viridis_c(option = "D") &
        theme(plot.title = element_text(size = 11))
      ggsave(file.path(OUT_DIR,
        paste0("tSNE_module_scores_page", page, ".pdf")),
        p_features, width = 15, height = 4 * ceiling(length(feat_subset) / 3))
    }
    cat("    Saved: tSNE_module_scores (", ceiling(n_row / 3), "pages)\n")
  }

  # --- 5f. Save final object ---
  saveRDS(sce, file.path(OUT_DIR, "sce_GSE192483_final.rds"))
  cat("  [5f] Saved: sce_GSE192483_final.rds\n")

  cat("  === Part 5 complete ===\n")
} else {
  cat("  [5] Skipping: annotated object not found. Run Part 2 first.\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n")
cat("============================================\n")
cat("  GSE192483 Analysis Complete\n")
cat("============================================\n")
cat("\n")

output_files <- list.files(OUT_DIR, pattern = "\\.(rds|pdf|csv)$")
if (length(output_files) > 0) {
  cat("Output files in", OUT_DIR, ":\n")
  for (f in sort(output_files)) {
    fpath <- file.path(normalizePath(OUT_DIR), f)
    fsize <- file.info(fpath)$size
    if (!is.na(fsize)) {
      if (fsize > 1e6) {
        fsize_str <- sprintf("%.1f MB", fsize / 1e6)
      } else if (fsize > 1e3) {
        fsize_str <- sprintf("%.1f KB", fsize / 1e3)
      } else {
        fsize_str <- sprintf("%d B", fsize)
      }
      cat(sprintf("  %-50s %s\n", f, fsize_str))
    } else {
      cat(sprintf("  %-50s\n", f))
    }
  }
} else {
  cat("  No output files yet. Run parts sequentially.\n")
}

cat("\n")
cat("Session info:\n")
sessionInfo()
