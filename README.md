# TB-TIGIT-NK Manuscript Code

Code repository for: **"Restoring NK Cell Effector Function through TIGIT Blockade Strengthens Host Defense against Mycobacterium tuberculosis"**

## Repository Structure

```
├── scRNA-seq/                          # Single-cell RNA-seq analysis
│   ├── 192483.R                        # GSE192483 lung scRNA-seq (QC, annotation, DEG)
│   ├── NK_full_analysis.R              # 1M-scBloodNL NK: TIGIT expression, effector DotPlot, module scores
│   ├── NK_metabolic_final.R            # NK metabolic gene timecourse (glycolysis/lactate/OXPHOS)
│   └── exhaustion_heatmap.R            # Parameterized TIGIT Low/High metabolic pathway heatmaps (4 sets)
├── bulkRNA-seq/                        # Bulk RNA-seq analysis
│   ├── α-TIGIT-修改.R                  # DEG (limma-voom), GSEA
│   ├── GSVA_1.R                        # GSVA pathway enrichment (KEGG + GO)
│   ├── WGCNA.R                         # Weighted gene co-expression network analysis
│   ├── ROSE.R                          # ROSE super-enhancer calling from H3K27ac ChIP-seq
│   └── AG分析SE.py                     # AlphaGenome deep-learning super-enhancer prediction
├── metabolomics/                       # Metabolomics analysis
│   └── TIGIT代谢组学测序分析.R          # DE, KEGG MSEA, Pathview, metabolite-gene correlation, DIABLO multi-omics
└── utils/                              # Shared utility functions
    ├── scRNA-seq.R                     # compute_pw_global_z(), process_tigit()
    └── go_kegg.R                       # run_go_enrichment(), run_kegg_enrichment(), run_gsea(), save_enrichment()
```

## Figure → Code Mapping

### Main Figures

| Figure Panel | Description | Code File |
|---|---|---|
| **Fig 1A** | Lung scRNA-seq UMAP, 13 cell types, tSNE | `scRNA-seq/192483.R` Parts 1–3 |
| **Fig 1B** | Donor-level Cyto/Exhaustion module scores (HC/FDR-LOW/FDR-HIGH) | `scRNA-seq/NK_full_analysis.R` Part 8 |
| **Fig 1C** | TIGIT expression in NK cells across UT/3h/24h | `scRNA-seq/NK_full_analysis.R` Part 1 |
| **Fig 3B** | GSVA GO heatmap (adj.P.Val < 0.05, |LogFC| > 1.0) | `bulkRNA-seq/GSVA_1.R` Parts 3–4 |
| **Fig 3C** | DEG volcano (cell migration & metabolism) | `bulkRNA-seq/α-TIGIT-修改.R` Part 2 |
| **Fig 3D** | GSEA (PI3K-AKT, mTOR, chemotaxis, cytotoxicity) | `bulkRNA-seq/α-TIGIT-修改.R` Part 3 |
| **Fig 5A** | Global z-score heatmap: mitochondrial energy + redox (TIGIT Low vs High, 3h + 24h) | `scRNA-seq/exhaustion_heatmap.R` (set: `mito_energy_redox`) |
| **Fig 5B** | GSEA: OXPHOS, mitochondrial pathways | `bulkRNA-seq/α-TIGIT-修改.R` Part 3 |
| **Fig 6A** | Volcano plot: 28 differential metabolites (VIP > 1, P < 0.05) | `metabolomics/TIGIT代谢组学测序分析.R` Part 3 |
| **Fig 6B** | KEGG enrichment of differential metabolites (MSEA, FDR < 0.05) | `metabolomics/TIGIT代谢组学测序分析.R` Part 6 |
| **Fig 6C** | scRNA-seq z-score heatmap: metabolomics-validated pathways (TIGIT Low vs High, 3h + 24h) | `scRNA-seq/exhaustion_heatmap.R` (set: `metabolomics`) |
| **Fig 6D** | Boxplots of 7 key metabolites + Pathview KEGG pathway maps | `metabolomics/TIGIT代谢组学测序分析.R` Parts 4, 8 |
| **Fig 6E** | Multi-omics integration (DIABLO: RNA + Metabolite blocks) | `metabolomics/TIGIT代谢组学测序分析.R` Part 10 |
| **Fig 6F** | Glycolysis/lactate/OXPHOS metabolic gene timecourse | `scRNA-seq/NK_metabolic_final.R` |
| **Fig 7A** | ROSE super-enhancer ranking ("hockey stick" plot) | `bulkRNA-seq/ROSE.R` |
| **Fig 7B** | AlphaGenome in silico H3K27ac prediction at TBX21 locus | `bulkRNA-seq/AG分析SE.py` |
| **Fig 8A** | Glycolytic/metabolic gene heatmap (1M-scBloodNL NK, UT/3h/24h) | `scRNA-seq/NK_metabolic_final.R` |

### Supplementary Figures

| Supp Panel | Description | Code File |
|---|---|---|
| **Supp 3A** | DotPlot: effector genes in TIGIT Low vs High NK (n=800/group) | `scRNA-seq/NK_full_analysis.R` Parts 3–4 |
| **Supp 4** | WGCNA: module-trait heatmap, MM-GS scatter, GO/KEGG enrichment, GSEA | `bulkRNA-seq/WGCNA.R` |

> **Note**: Figures 2, 4, 5C–Q, 6F–N, 7C–O, 8B–P, 9, Supp 1–2, Supp 3B–C, Supp 5–8 contain wet-lab experimental data (flow cytometry, qPCR, western blot, confocal microscopy, CFU assays) and do not have corresponding computational code in this repository.

## Data Dependencies

All scripts expect input data in a `./data/` subdirectory and write outputs to `./output/`.

| Dataset | Accession / Source |
|---|---|
| 1M-scBloodNL MTB PBMC scRNA-seq | EGA: EGAS00001005376 |
| GSE192483 (Lung TB scRNA-seq) | GEO |
| GSE227691 (Healthy lung scRNA-seq) | GEO |
| ENCFF119TOL (Primary NK H3K27ac ChIP-seq) | ENCODE |
| 114 metabolic pathway gene sets | `114条代谢通路基因.xlsx` (included in data directory) |
| Bulk RNA-seq (NK-92 anti-TIGIT vs IgG) | Contact corresponding author |
| Targeted metabolomics (NK-92 anti-TIGIT vs IgG) | Contact corresponding author |

## Usage

Expensive one-time computations (QC clustering, DE, WGCNA network construction) are wrapped in `if (FALSE)` blocks — set to `TRUE` on first run.

### R scripts
```r
source("scRNA-seq/NK_full_analysis.R")
source("bulkRNA-seq/WGCNA.R")
source("metabolomics/TIGIT代谢组学测序分析.R")
```

### Python
```bash
python bulkRNA-seq/AG分析SE.py
```

## Dependencies

- **R**: Seurat v5, harmony, limma, edgeR, DESeq2, WGCNA, GSVA, clusterProfiler, enrichplot, pathview, pheatmap, ggplot2, mixOmics, ropls, randomForest, corrplot, psych, ConsensusClusterPlus
- **Python**: alphagenome, numpy, pandas, matplotlib
- **Bioconductor**: org.Hs.eg.db, msigdbr, KEGGREST
