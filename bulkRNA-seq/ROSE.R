# ======================================================================
# ROSE.R — ROSE (Rank Ordering of Super-Enhancers)
# H3K27ac ChIP-seq peak stitching + bigWig quantification
# ======================================================================

# ---- Setup -----------------------------------------------------------

base_dir <- getwd()
DATA_DIR <- base_dir  # Change this to your data directory if needed

cat("===== ROSE Super-Enhancer Calling =====\n")
cat("Base directory:", base_dir, "\n")
cat("Data directory:", DATA_DIR, "\n\n")

dir.create(file.path(base_dir, "output"), showWarnings = FALSE)

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ---- Input files -----------------------------------------------------

peak_file <- file.path(DATA_DIR,
  "GSM6727286_Ishikawa_8hrDMSO_H3K27Ac_rep1_stringent_peaks.narrowPeak.gz")
bw_file   <- file.path(DATA_DIR,
  "GSM6727286_Ishikawa_8hrDMSO_H3K27Ac_rep1_stringent.bw")
tss_file  <- file.path(DATA_DIR, "hg38_refseq_tss.txt")

# ======================================================================
# STEP 1: Load peaks
# ======================================================================

cat(">>> Step 1: Load peaks\n")

peaks <- fread(
  peak_file, header = TRUE,
  col.names = c("chr", "start", "end", "name", "score", "strand",
                "signalValue", "pValue", "qValue", "peak")
)

peaks_gr <- GRanges(peaks$chr, IRanges(peaks$start, peaks$end))
cat("Total peaks loaded:", nrow(peaks), "\n")

# ======================================================================
# STEP 2: Load unique TSS (from refGene)
# ======================================================================

cat(">>> Step 2: Load TSS (from hg38_refseq_tss.txt)\n")

tss <- fread(tss_file, header = FALSE,
             col.names = c("chr", "tss_start", "tss_end", "gene"))

# One TSS per gene: choose median TSS
tss_unique <- tss[, .(tss = median(tss_start)), by = .(chr, gene)]

tss_gr <- GRanges(
  seqnames = tss_unique$chr,
  ranges   = IRanges(tss_unique$tss - 2500,
                     tss_unique$tss + 2500)
)

cat("Unique TSS:", nrow(tss_unique), "\n")

# ======================================================================
# STEP 3: Remove promoter-proximal peaks (+/- 2.5 kb around TSS)
# ======================================================================

cat(">>> Step 3: Remove promoter peaks\n")

is_promoter <- overlapsAny(peaks_gr, tss_gr)
enh <- peaks[!is_promoter]
enh_gr <- GRanges(enh$chr, IRanges(enh$start, enh$end))

cat("Remaining enhancer peaks:", nrow(enh), "\n")

# ======================================================================
# STEP 4: Quantify enhancer signal from bigWig
# ======================================================================

cat(">>> Step 4: Quantify enhancer signal from bigWig\n")

bw <- import(bw_file)

ov <- findOverlaps(enh_gr, bw)

dt <- data.table(
  idx = queryHits(ov),
  cov = width(pintersect(enh_gr[queryHits(ov)],
                         bw[subjectHits(ov)])) *
    bw$score[subjectHits(ov)]
)

signal_dt <- dt[, .(signal = sum(cov)), by = idx]

enh$signal <- 0
enh$signal[signal_dt$idx] <- signal_dt$signal

# ======================================================================
# STEP 5: Stitch enhancers within 12.5 kb (ROSE standard)
# ======================================================================

cat(">>> Step 5: Stitching enhancers\n")

enh <- enh[order(chr, start)]

stitched <- list()
cur_chr <- enh$chr[1]
cur_start <- enh$start[1]
cur_end <- enh$end[1]
cur_signal <- enh$signal[1]

for (i in 2:nrow(enh)) {
  if (enh$chr[i] == cur_chr &&
      (enh$start[i] - cur_end) <= 12500) {
    cur_end <- max(cur_end, enh$end[i])
    cur_signal <- cur_signal + enh$signal[i]
  } else {
    stitched[[length(stitched) + 1]] <- data.table(
      chr = cur_chr,
      stitched_start = cur_start,
      stitched_end = cur_end,
      signal = cur_signal
    )
    cur_chr <- enh$chr[i]
    cur_start <- enh$start[i]
    cur_end <- enh$end[i]
    cur_signal <- enh$signal[i]
  }
}

stitched[[length(stitched) + 1]] <- data.table(
  chr = cur_chr,
  stitched_start = cur_start,
  stitched_end = cur_end,
  signal = cur_signal
)

stitched <- rbindlist(stitched)
cat("Stitched enhancer loci:", nrow(stitched), "\n")

# ======================================================================
# STEP 6: Super-enhancer ranking (ROSE tangent-line method)
# ======================================================================

cat(">>> Step 6: Ranking super-enhancers\n")

stitched <- stitched[order(signal)]
stitched[, rank := 1:.N]

stitched[, norm_signal := signal / max(signal)]
stitched[, norm_rank := rank / max(rank)]
stitched[, dist := abs(norm_signal - norm_rank)]
cutoff_idx <- stitched[which.max(dist), rank]
super_enhancers <- stitched[rank >= cutoff_idx]

cat("Super-enhancers called:", nrow(super_enhancers), "\n")

# ======================================================================
# STEP 7: Assign nearest gene to each super-enhancer
# ======================================================================

cat(">>> Step 7: Assign nearest gene\n")

se_gr <- GRanges(
  super_enhancers$chr,
  IRanges(super_enhancers$stitched_start,
          super_enhancers$stitched_end)
)

hits <- distanceToNearest(se_gr, tss_gr)

super_enhancers$gene <- tss_unique$gene[subjectHits(hits)]
super_enhancers$dist <- mcols(hits)$distance

# ======================================================================
# STEP 8: Plot super-enhancer ranking curve
# ======================================================================

cat(">>> Step 8: Plotting\n")

genes_to_label <- c("TBX21", "EOMES", "TOX", "PRDM1", "STAT4", "STAT5B",
                    "RUNX3", "BATF", "IKZF3", "GZMB", "PRF1", "GZMK",
                    "FASLG", "CTLA4", "CD226", "IFNG")

p <- ggplot(stitched, aes(rank, signal)) +
  geom_line(size = 1.2, color = "grey40") +
  geom_point(data = super_enhancers, color = "red", size = 2) +
  geom_vline(xintercept = cutoff_idx, linetype = "dashed") +
  theme_classic(base_size = 16) +
  geom_text_repel(
    data = super_enhancers[gene %in% genes_to_label],
    aes(label = gene),
    fontface = "bold",
    size = 5,
    force = 25,
    max.overlaps = Inf
  ) +
  labs(
    title = "Super-Enhancer Landscape (ROSE, bigWig-based)",
    x = "Enhancer Rank",
    y = "H3K27ac Signal (bigWig x width)"
  )

print(p)

ggsave(
  plot = p,
  filename = file.path(base_dir, "output", "NK_SuperEnhancers_ROSE_plot.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# ======================================================================
# STEP 9: Save results
# ======================================================================

cat(">>> Step 9: Saving results\n")

write.csv(super_enhancers,
          file.path(base_dir, "output", "NK_SuperEnhancers_ROSE.csv"),
          row.names = FALSE)

# Human-readable sorted version
super_enhancers_pretty <- copy(super_enhancers)
super_enhancers_pretty <- super_enhancers_pretty[order(-signal)]
super_enhancers_pretty[, SE_rank := 1:.N]

write.csv(super_enhancers_pretty,
          file.path(base_dir, "output", "SuperEnhancers_final_human_readable.csv"),
          row.names = FALSE)

# ======================================================================
# Summary
# ======================================================================

cat("\n===== ROSE complete =====\n")
cat("Super-enhancers called:", nrow(super_enhancers), "\n")
cat("Output directory:", file.path(base_dir, "output"), "\n\n")

# Print genes of interest
print(super_enhancers[grepl("TOX|EOMES|ZEB2|IKZF3", gene)])
head(super_enhancers)

cat("\nAll outputs saved to: ", file.path(base_dir, "output"), "\n")
