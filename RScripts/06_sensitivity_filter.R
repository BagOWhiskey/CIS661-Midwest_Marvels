# =============================================================================
# 06_sensitivity_filter.R
# Proposal section F: how much does the DE count depend on OUR probe-filtering
# choice? (The paper does not say how it filtered.) Reruns the Aim 1 limma
# analysis under several filter settings and tabulates the results.
#
# Input:  data/processed/eset_rma_annotated.rds (from 02_qc_normalize.R)
# Output: results/sensitivity_filter.csv
# Runtime: about a minute; no plots, so memory use stays low.
# =============================================================================

# Work from the project folder even if RStudio has the parent project open
if (!file.exists(file.path("R", "_config.R")) &&
    file.exists(file.path("mir375_replication", "R", "_config.R"))) setwd("mir375_replication")

source(file.path("R", "_config.R"))
suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
})
set.seed(SEED)

eset_all <- readRDS(ESET_FILE)
eset_all$group <- factor(eset_all$group, levels = GROUP_LEVELS)

# Filter settings to compare. The row matching _config.R is our main analysis.
settings <- data.frame(
  setting           = c("No filtering", "Annotated only",
                        "Annotated + low-expr 10%", "Annotated + low-expr 20% (main)",
                        "Annotated + low-expr 30%"),
  annotated_only    = c(FALSE, TRUE, TRUE, TRUE, TRUE),
  low_expr_quantile = c(NA,    NA,   0.10, 0.20, 0.30)
)

run_setting <- function(annotated_only, low_expr_quantile) {
  eset <- eset_all
  if (annotated_only) eset <- eset[!is.na(fData(eset)$SYMBOL), ]
  if (!is.na(low_expr_quantile)) {
    cutoff <- quantile(exprs(eset), low_expr_quantile)
    eset   <- eset[rowSums(exprs(eset) > cutoff) >= LOW_EXPR_MIN_ARRAYS, ]
  }

  design <- model.matrix(~ 0 + group, data = pData(eset))
  colnames(design) <- levels(eset$group)
  contrasts <- makeContrasts(Dox_0 - WT, Dox_10 - WT, Dox_30 - WT,
                             Dox_10 - Dox_0, Dox_30 - Dox_0, Dox_30 - Dox_10,
                             levels = design)
  fit2 <- eBayes(contrasts.fit(lmFit(eset, design), contrasts))
  dt   <- decideTests(fit2, method = "separate", adjust.method = "BH", p.value = ALPHA)

  # The six highlighted genes, 30 vs 0 ng/mL (column 5 = Dox_30 - Dox_0)
  tt   <- topTable(fit2, coef = 5, number = Inf, adjust.method = "BH")
  hit  <- tt$SYMBOL %in% HIGHLIGHT_GENES & tt$logFC < 0 & tt$adj.P.Val < ALPHA
  n_hi <- length(unique(tt$SYMBOL[hit]))

  data.frame(
    probesets_tested    = nrow(eset),
    de_any_contrast     = sum(rowSums(dt != 0) > 0),
    de_30_vs_0          = sum(dt[, 5] != 0),
    highlighted_down_30 = paste0(n_hi, " of 6"),
    pct_of_paper_7656   = round(100 * sum(rowSums(dt != 0) > 0) / PAPER$n_de_any_contrast)
  )
}

results <- do.call(rbind, lapply(seq_len(nrow(settings)), function(i) {
  message("Running: ", settings$setting[i])
  cbind(settings[i, ], run_setting(settings$annotated_only[i], settings$low_expr_quantile[i]))
}))
rownames(results) <- NULL

print(results[, c("setting", "probesets_tested", "de_any_contrast",
                  "de_30_vs_0", "highlighted_down_30", "pct_of_paper_7656")], row.names = FALSE)
write.csv(results, file.path(DIR_RESULTS, "sensitivity_filter.csv"), row.names = FALSE)
message("Saved results/sensitivity_filter.csv (paper reports ", PAPER$n_de_any_contrast,
        " DE transcripts)")
