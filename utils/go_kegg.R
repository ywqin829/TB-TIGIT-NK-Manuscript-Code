# ==================== utils/go_kegg.R ====================
# Shared GO/KEGG enrichment wrapper functions
# TB-TIGIT-NK manuscript
# ==========================================================

#' Run GO enrichment analysis with clusterProfiler
#'
#' @param gene_list Character vector of gene SYMBOLs
#' @param ont       Ontology: "BP", "MF", "CC", or "all"
#' @param OrgDb     org.*.eg.db package name
#' @param keyType   Key type (default "SYMBOL")
#' @return enrichResult object
run_go_enrichment <- function(gene_list, ont = "all",
                               OrgDb = "org.Hs.eg.db",
                               keyType = "SYMBOL") {
  ids <- tryCatch(
    bitr(gene_list, fromType = keyType, toType = "ENTREZID", OrgDb = OrgDb),
    error = function(e) {
      cat("  WARNING: bitr conversion failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(ids) || nrow(ids) == 0) {
    cat("  No valid ENTREZ IDs found.\n")
    return(NULL)
  }
  ego <- enrichGO(
    gene = ids$ENTREZID,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  ego
}

#' Run KEGG enrichment analysis with clusterProfiler
#'
#' @param gene_list Character vector of gene SYMBOLs
#' @param organism  KEGG organism code (default "hsa")
#' @param OrgDb     org.*.eg.db package name
#' @return enrichResult object
run_kegg_enrichment <- function(gene_list, organism = "hsa",
                                 OrgDb = "org.Hs.eg.db") {
  ids <- tryCatch(
    bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb),
    error = function(e) {
      cat("  WARNING: bitr conversion failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(ids) || nrow(ids) == 0) {
    cat("  No valid ENTREZ IDs found.\n")
    return(NULL)
  }
  kk <- enrichKEGG(
    gene = ids$ENTREZID,
    organism = organism,
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  kk
}

#' Run GSEA (GO and KEGG)
#'
#' @param geneList Named numeric vector (ENTREZID = logFC), sorted decreasing
#' @param OrgDb    org.*.eg.db package name
#' @return List with $go and $kegg results
run_gsea <- function(geneList, OrgDb = "org.Hs.eg.db") {
  go_result <- tryCatch(
    gseGO(geneList, OrgDb = OrgDb, keyType = "ENTREZID",
          ont = "BP", eps = 0, pvalueCutoff = 0.05, seed = 42),
    error = function(e) {
      cat("  GSEA GO failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  kegg_result <- tryCatch(
    gseKEGG(geneList, pvalueCutoff = 0.05, eps = 0, seed = 42),
    error = function(e) {
      cat("  GSEA KEGG failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  list(go = go_result, kegg = kegg_result)
}

#' Save enrichment results to CSV
#'
#' @param ego enrichResult object
#' @param file_path Output CSV path
save_enrichment <- function(ego, file_path) {
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file_path, row.names = FALSE)
    cat("  Saved:", file_path, "\n")
  } else {
    cat("  No results to save for:", file_path, "\n")
  }
}
