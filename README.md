# TB-TIGIT-NK Manuscript Code

Code for single-cell and bulk transcriptomic analysis of TIGIT+ NK cell metabolic exhaustion in tuberculosis.

## Overview

This repository contains analysis scripts for the TB-TIGIT-NK manuscript, which investigates the metabolic and functional characteristics of TIGIT-expressing NK cells in tuberculosis patients. The analysis integrates scRNA-seq, bulk RNA-seq, metabolomics, and super-enhancer profiling.

## Repository Structure

```
├── scRNA-seq/               # Single-cell RNA-seq analysis
│   ├── 192483.R             # GSE192483 pre-processing + T-cell DEG
│   ├── NK_full_analysis.R   # NK cell main analysis (TIGIT expression, pathway heatmaps, cytotoxicity/exhaustion modules)
│   ├── NK_metabolic_final.R # Gene-level metabolic enzyme/transporter analysis
│   └── exhaustion_heatmap.R # Parameterized TIGIT Low/High pathway heatmaps (4 pathway sets)
├── bulkRNA-seq/              # Bulk RNA-seq analysis
│   ├── α-TIGIT-修改.R        # Differential expression (limma-voom), GO/KEGG, GSEA, CIBERSORT
│   ├── GSVA_1.R             # GSVA pathway enrichment (KEGG + GO)
│   ├── WGCNA.R              # Weighted gene co-expression network analysis
│   ├── ROSE.R               # ROSE super-enhancer calling from ChIP-seq
│   └── AG分析SE.py          # AlphaGenome-based super-enhancer prediction
├── metabolomics/             # Metabolomics analysis
│   └── TIGIT代谢组学测序分析.R  # Metabolite differential expression, KEGG enrichment, CIBERSORT
└── utils/                    # Shared utility functions
    ├── scRNA-seq.R           # compute_pw_global_z(), process_tigit(), build_tigit_heatmap()
    └── go_kegg.R             # run_go_enrichment(), run_kegg_enrichment(), run_gsea(), save_enrichment()
```

## Data Dependencies

- **MTB timecourse scRNA-seq**: ~358K PBMC from 120 donors (in-house dataset, contact corresponding author for access)
- **GSE192483**: Lung scRNA-seq from TB patients (GEO)
- **GSE40553**: Bulk microarray data (GEO)
- **Pathway gene sets**: `114条代谢通路基因.xlsx` — 114 metabolic pathway gene list

## Usage

All scripts expect data files in a `./data/` subdirectory and output figures/CSVs to `./output/`. Expensive one-time computations (QC clustering, differential expression, WGCNA network construction) are wrapped in `if (FALSE)` blocks — set to `TRUE` on first run, then `FALSE` for subsequent re-runs using cached results.

### R scripts
```r
# Example: NK cell main analysis
source("scRNA-seq/NK_full_analysis.R")
```

### Python script
```bash
python bulkRNA-seq/AG分析SE.py
```

## Key Analyses

- **TIGIT+ NK cell metabolic exhaustion**: Pathway-level z-score analysis comparing TIGIT High vs Low NK cells across 3h and 24h MTB stimulation
- **Super-enhancer prediction**: ROSE (R) + AlphaGenome (Python) for TBX21 super-enhancer characterization
- **Immune deconvolution**: CIBERSORT-based estimation of immune cell composition
- **Co-expression networks**: WGCNA for module-trait correlation analysis

## Output

All figures and tables are saved to the `output/` directory as PDFs and CSVs.

## Session Info

Analysis performed with R 4.x, Seurat v5, Bioconductor packages, and Python 3.x.

## Related Repository

- [TF-TB-Manuscript-Code](https://github.com/ywqin829/TF-TB-Manuscript-Code) — TBX21-PAX5 axis in TB resolution (*Respiratory Medicine*, 2025)
