# =============================================================================
# CIS 661 Group Project
# Group 1: Joe Baczkowski, Minh Tran, Faizaan Khan, Patrick Schneider
# 02_qc_normalize.R
# Proposal section E1: reads the 16 Clariom D CEL files, makes raw-data QC
# plots, normalizes with RMA at the gene (transcript-cluster) level, annotates
# probesets with gene symbols, and makes post-normalization QC plots.
#
# Output: data/processed/eset_rma_annotated.rds (normalized ExpressionSet)
# Memory: reading 16 Clariom D arrays needs several GB of RAM; use the GVSU
#         RStudio Server rather than a small laptop if it runs out of memory.
# =============================================================================

# Work from the project folder even if RStudio has the parent project open
if (!file.exists(file.path("R", "_config.R")) &&
    file.exists(file.path("mir375_replication", "R", "_config.R"))) setwd("mir375_replication")

source(file.path("R", "_config.R"))
suppressPackageStartupMessages({
  if (LOW_MEMORY_MODE) library(ff)   # must be loaded before oligo to take effect
  library(oligo)
  library(Biobase)
  library(AnnotationDbi)
  library(clariomdhumantranscriptcluster.db)
  library(ggplot2)
})
set.seed(SEED)

# Low-memory mode: oligo stores the probe-level intensities in files under FF_DIR
# and processes a few arrays at a time, which keeps the session under 4 GB
if (LOW_MEMORY_MODE) {
  dir.create(FF_DIR, recursive = TRUE, showWarnings = FALSE)
  ldPath(FF_DIR)
  ocSamples(4)        # arrays processed per batch
  ocProbesets(50000)  # probesets processed per batch
  message("Low-memory mode on (ff files in ", FF_DIR, ")")
}

sample_sheet <- read.csv(SAMPLE_SHEET_FILE, stringsAsFactors = FALSE)
sample_sheet$group <- factor(sample_sheet$group, levels = GROUP_LEVELS)
group_cols <- GROUP_COLORS[as.character(sample_sheet$group)]

save_png <- function(file, expr, width = 1600, height = 1000) {
  png(file.path(DIR_FIGURES, file), width = width, height = height, res = 150)
  on.exit(dev.off())
  expr
}

# ---- 1. Read raw CEL files ------------------------------------------------------------
raw <- read.celfiles(sample_sheet$cel_file)        # uses pd.clariom.d.human
sampleNames(raw) <- sample_sheet$sample_id

# Scan dates, to check for batch effects in the PCA below
# (CEL headers store e.g. "2024-03-01T10:22:11Z"; keep only the date part)
scan_dates <- sub("[T ].*$", "", as.character(protocolData(raw)$dates))
sample_sheet$scan_date <- scan_dates
message("Scan dates:"); print(table(sample_sheet$group, sample_sheet$scan_date))

# ---- 2. Raw-data QC -------------------------------------------------------------------------
# (A failed plot should not stop the normalization, so each one is wrapped)
tryCatch(save_png("qc_01_raw_boxplot.png",
                  oligo::boxplot(raw, target = "core", las = 2, col = group_cols,
                                 main = "Raw log2 intensities (before RMA)")),
         error = function(e) warning("Raw boxplot skipped: ", conditionMessage(e)))
tryCatch(save_png("qc_02_raw_density.png",
                  oligo::hist(raw, target = "core", col = group_cols,
                              main = "Raw intensity distributions")),
         error = function(e) warning("Raw density plot skipped: ", conditionMessage(e)))

# ---- 3. RMA normalization (background correction + quantile normalization +
#         median-polish summarization) at the transcript-cluster ("gene") level ----------
eset <- rma(raw, target = "core")
rm(raw); invisible(gc())

# The gene-level result is small: make sure it is an ordinary in-memory matrix,
# then delete the temporary probe-level ff files
if (!is.matrix(exprs(eset))) exprs(eset) <- as.matrix(exprs(eset)[, ])
if (LOW_MEMORY_MODE) unlink(FF_DIR, recursive = TRUE)
message("RMA done: ", nrow(eset), " probesets x ", ncol(eset), " arrays")

# ---- 4. Gene annotation -------------------------------------------------------------------
ann <- AnnotationDbi::select(clariomdhumantranscriptcluster.db,
                             keys    = featureNames(eset),
                             columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                             keytype = "PROBEID")
# A few probesets map to several genes: keep the first mapping, but record how many
n_multi <- sum(duplicated(ann$PROBEID))
ann <- ann[!duplicated(ann$PROBEID), ]
ann <- ann[match(featureNames(eset), ann$PROBEID), ]
rownames(ann) <- featureNames(eset)
fData(eset) <- ann
message(sum(!is.na(ann$SYMBOL)), " of ", nrow(ann),
        " probesets have a gene symbol (", n_multi, " extra multi-gene mappings dropped)")

# Sample information travels with the data from here on
pd <- sample_sheet
rownames(pd) <- sample_sheet$sample_id
pData(eset) <- pd[sampleNames(eset), ]

# ---- 5. Post-normalization QC ----------------------------------------------------------------
expr <- exprs(eset)

save_png("qc_03_rma_boxplot.png", {
  par(mar = c(7, 4, 3, 1))
  boxplot(expr, las = 2, col = group_cols, outline = FALSE,
          ylab = "log2 expression", main = "After RMA normalization")
})

# PCA on annotated probesets (all samples should cluster by group, not by scan date)
keep <- !is.na(fData(eset)$SYMBOL)
pca  <- prcomp(t(expr[keep, ]), center = TRUE, scale. = FALSE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
pca_df <- data.frame(pca$x[, 1:2], pData(eset))
p <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = factor(scan_date))) +
  geom_point(size = 3.5) +
  scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS) +
  labs(x = paste0("PC1 (", var_pct[1], "%)"), y = paste0("PC2 (", var_pct[2], "%)"),
       color = "Group", shape = "Scan date",
       title = "PCA of RMA-normalized arrays (GSE272451)") +
  theme_bw()
ggsave(file.path(DIR_FIGURES, "qc_04_pca.png"), p, width = 7, height = 5, dpi = 150)

# Sample-to-sample correlation (an outlier array shows up as a pale row/column)
cor_mat <- cor(expr[keep, ])
save_png("qc_05_sample_correlation.png",
         heatmap(cor_mat, symm = TRUE, ColSideColors = group_cols,
                 RowSideColors = group_cols, main = "Sample correlation"),
         width = 1300, height = 1300)
write.csv(round(cor_mat, 4), file.path(DIR_RESULTS, "qc_sample_correlation.csv"))

# ---- 6. Save -----------------------------------------------------------------------------------
saveRDS(eset, ESET_FILE)
message("Saved normalized ExpressionSet to ", ESET_FILE)
writeLines(capture.output(sessionInfo()), file.path(DIR_RESULTS, "sessionInfo_02.txt"))
