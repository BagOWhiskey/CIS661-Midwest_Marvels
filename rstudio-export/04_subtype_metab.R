# =============================================================================
# 04_subtype_metab.R
# Aim 2: does miR-375 shift C4-2B cells from MetB-like to MetA-like? (paper Fig. 3C)
#
# The paper used CIBERSORTx with reference profiles built from 157 MetA-C transcripts
# (Thysell et al., 2022). Those reference profiles were never published, so a true
# CIBERSORTx rerun is not possible. This script therefore:
#   1. Writes the mixture file CIBERSORTx needs, in case the reference can be obtained
#      (e.g. from the authors), so that the original method can still be run on the web.
#   2. Uses the proposal's fallback: scores every array for MetA, MetB and MetC with the
#      mean z-score of each subtype's genes from Table S2, then tests whether the
#      MetA - MetB score rises with miR-375 dose (unpaired t-tests, as in Fig. 3C).
#
# Input:  data/processed/eset_rma_annotated.rds (from 02), Table S2 (see _config.R)
# Output: results/cibersortx_mixture_GSE272451.txt, results/aim2_*.csv,
#         figures/aim2_subtype_scores.png, figures/aim2_signature_heatmap.png
# NOTE:   the scores are RELATIVE to these 16 arrays (z-scores). They show direction and
#         significance of a shift, not the paper's "% similarity", which needs the
#         CIBERSORTx reference. Report this as a deviation from the original method.
# =============================================================================

# Work from the project folder even if RStudio has the parent project open
if (!file.exists(file.path("R", "_config.R")) &&
    file.exists(file.path("mir375_replication", "R", "_config.R"))) setwd("mir375_replication")

source(file.path("R", "_config.R"))
suppressPackageStartupMessages({
  library(Biobase)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})
set.seed(SEED)

eset <- readRDS(ESET_FILE)
eset$group <- factor(eset$group, levels = GROUP_LEVELS)

# ---- 1. One row per gene (log2), and the CIBERSORTx mixture file ----------------------------
# Several probesets can map to one gene: keep the one with the highest mean expression
annotated <- !is.na(fData(eset)$SYMBOL)
expr_ann  <- exprs(eset)[annotated, ]
sym_ann   <- fData(eset)$SYMBOL[annotated]
best      <- order(sym_ann, -rowMeans(expr_ann))
best      <- best[!duplicated(sym_ann[best])]
gene_expr <- expr_ann[best, ]
rownames(gene_expr) <- sym_ann[best]
message(nrow(gene_expr), " genes after collapsing probesets to gene symbols")

# CIBERSORTx expects non-log expression, genes in rows, first column = gene symbol
mixture <- data.frame(GeneSymbol = rownames(gene_expr), 2^gene_expr, check.names = FALSE)
write.table(mixture, file.path(DIR_RESULTS, "cibersortx_mixture_GSE272451.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("CIBERSORTx mixture file written to results/cibersortx_mixture_GSE272451.txt")

# ---- 2. Read Table S2 (157 MetA-C transcripts) -----------------------------------------------
if (!file.exists(SIGNATURE_FILE)) {
  stop("Signature file not found: ", SIGNATURE_FILE, "\n",
       "Download Table S2 from Thysell et al. (2022), https://doi.org/10.1002/1878-0261.13158 ",
       "(Supporting Information), save it under that name, and rerun this script.")
}

# Read everything as text without assuming where the header is (supplementary tables
# often have a title line above the real column names)
read_any <- function(path) {
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
    raw <- as.data.frame(suppressMessages(
      readxl::read_excel(path, col_names = FALSE, col_types = "text", .name_repair = "minimal")))
  } else {
    raw <- read.csv(path, header = FALSE, colClasses = "character", check.names = FALSE)
  }
  raw[is.na(raw)] <- ""
  raw
}
to_subtype <- function(x) {
  x <- toupper(gsub("[^A-Za-z]", "", x))
  x <- sub("^(METASTASISSUBTYPE|SUBTYPE|METASTASIS|MET)", "", x)
  ifelse(x %in% c("A", "B", "C"), paste0("Met", x), NA_character_)
}

raw <- read_any(SIGNATURE_FILE)
first_rows <- raw[seq_len(min(15, nrow(raw))), , drop = FALSE]

# Two layouts are handled:
#  (a) "wide": one column per subtype, headed MetA | MetB | MetC, genes listed below
#      (this is how the published Table S2 is laid out)
#  (b) "long": a gene-symbol column plus a subtype column
n_subtype_labels <- function(r) sum(!is.na(to_subtype(r)) & nchar(trimws(r)) <= 12)
wide_row <- which(apply(first_rows, 1, n_subtype_labels) >= 2)[1]

if (!is.na(wide_row) && is.null(SIGNATURE_GENE_COL)) {
  subtype_of_col <- to_subtype(unlist(raw[wide_row, ]))
  below <- raw[-seq_len(wide_row), , drop = FALSE]
  signature <- do.call(rbind, lapply(which(!is.na(subtype_of_col)), function(j)
    data.frame(symbol = trimws(below[[j]]), subtype = subtype_of_col[j], stringsAsFactors = FALSE)))
  message("Table S2 layout: one column per subtype (row ", wide_row, ")")
} else {
  # The header is the first row with at least two filled cells, one of them a short label
  # mentioning "gene" or "symbol" (a long title sentence above the table does not count)
  is_header <- function(r) {
    r <- trimws(r)
    sum(nzchar(r)) >= 2 && any(nchar(r) <= 30 & grepl("symbol|gene", r, ignore.case = TRUE))
  }
  header_row <- which(apply(first_rows, 1, is_header))[1]
  if (is.na(header_row)) stop("Could not find the table header in ", SIGNATURE_FILE,
                              " (expected MetA/MetB/MetC columns, or a gene-symbol column)")
  sig_tab <- raw[-seq_len(header_row), , drop = FALSE]
  names(sig_tab) <- make.unique(trimws(unlist(raw[header_row, ])))
  sig_tab <- sig_tab[apply(sig_tab, 1, function(r) any(nzchar(r))), , drop = FALSE]

  gene_col <- SIGNATURE_GENE_COL
  if (is.null(gene_col)) {
    gene_col <- c(grep("symbol", names(sig_tab), ignore.case = TRUE, value = TRUE),
                  grep("gene",   names(sig_tab), ignore.case = TRUE, value = TRUE))[1]
  }
  subtype_col <- SIGNATURE_SUBTYPE_COL
  if (is.null(subtype_col)) {
    share <- vapply(sig_tab, function(v) mean(!is.na(to_subtype(v))), numeric(1))
    subtype_col <- if (max(share) >= 0.8) names(sig_tab)[which.max(share)] else NA
  }
  if (is.na(gene_col) || is.na(subtype_col)) {
    message("Columns found in Table S2: ", paste(names(sig_tab), collapse = " | "))
    print(head(sig_tab))
    stop("Could not identify the gene-symbol and/or subtype column. Set SIGNATURE_GENE_COL ",
         "and SIGNATURE_SUBTYPE_COL in R/_config.R to the exact column names above.")
  }
  message("Table S2: gene column = '", gene_col, "', subtype column = '", subtype_col, "'")
  signature <- data.frame(symbol  = trimws(sig_tab[[gene_col]]),
                          subtype = to_subtype(sig_tab[[subtype_col]]),
                          stringsAsFactors = FALSE)
}
signature <- signature[nzchar(signature$symbol) & !is.na(signature$subtype), ]
signature <- signature[!duplicated(signature), ]
message(nrow(signature), " signature transcripts read (paper: 157); per subtype:")
print(table(signature$subtype))

# ---- 3. Match signature genes to the array (rescuing old gene names via aliases) --------------
signature$matched <- ifelse(signature$symbol %in% rownames(gene_expr), signature$symbol, NA)
unmatched <- is.na(signature$matched)
if (any(unmatched)) {
  # (mapIds stops with an error if none of the names is a known alias, hence tryCatch)
  alias <- tryCatch(unname(suppressMessages(mapIds(org.Hs.eg.db,
                      keys = signature$symbol[unmatched], column = "SYMBOL",
                      keytype = "ALIAS", multiVals = "first"))),
                    error = function(e) rep(NA_character_, sum(unmatched)))
  signature$matched[unmatched] <- ifelse(!is.na(alias) & alias %in% rownames(gene_expr), alias, NA)
}
signature$how <- ifelse(is.na(signature$matched), "not on array",
                        ifelse(signature$matched == signature$symbol, "symbol", "alias"))
write.csv(signature, file.path(DIR_RESULTS, "aim2_signature_genes_matched.csv"), row.names = FALSE)
used <- signature[!is.na(signature$matched) & !duplicated(signature$matched), ]
message(nrow(used), " of ", nrow(signature), " signature genes found on the array:")
print(table(used$subtype))
if (any(table(factor(used$subtype, levels = c("MetA", "MetB"))) < 5)) {
  stop("Fewer than 5 MetA or MetB genes matched; check the signature file.")
}

# ---- 4. Subtype scores: mean row z-score of each subtype's genes ---------------------------
z <- t(scale(t(gene_expr[used$matched, , drop = FALSE])))
subtypes <- sort(unique(used$subtype))
scores <- sapply(subtypes, function(s) colMeans(z[used$subtype == s, , drop = FALSE]))
scores <- data.frame(sample = colnames(z), group = eset$group, scores, check.names = FALSE)
scores$MetA_minus_MetB <- scores$MetA - scores$MetB
write.csv(scores, file.path(DIR_RESULTS, "aim2_subtype_scores_per_sample.csv"), row.names = FALSE)

score_cols <- c(subtypes, "MetA_minus_MetB")
group_summary <- do.call(rbind, lapply(score_cols, function(sc)
  data.frame(score = sc, group = GROUP_LEVELS,
             mean = tapply(scores[[sc]], scores$group, mean),
             sd   = tapply(scores[[sc]], scores$group, sd))))
rownames(group_summary) <- NULL
print(group_summary, digits = 2, row.names = FALSE)
write.csv(group_summary, file.path(DIR_RESULTS, "aim2_subtype_scores_by_group.csv"), row.names = FALSE)

# ---- 5. Unpaired t-tests (as in the paper's Fig. 3C) ----------------------------------------
comparisons <- list(c("Dox_0", "WT"), c("Dox_10", "Dox_0"), c("Dox_30", "Dox_0"), c("Dox_30", "WT"))
tests <- do.call(rbind, lapply(score_cols, function(sc) do.call(rbind, lapply(comparisons, function(cmp) {
  a <- scores[[sc]][scores$group == cmp[1]]
  b <- scores[[sc]][scores$group == cmp[2]]
  tt <- t.test(a, b, var.equal = TRUE)   # Student's unpaired t-test, n = 4 per group
  data.frame(score = sc, comparison = paste(cmp[1], "vs", cmp[2]),
             difference = mean(a) - mean(b), t = unname(tt$statistic), p_value = tt$p.value)
}))))
tests$significant <- tests$p_value < ALPHA
print(tests, digits = 3, row.names = FALSE)
write.csv(tests, file.path(DIR_RESULTS, "aim2_ttests.csv"), row.names = FALSE)

# Comparison with the paper (proposal section E4)
diff_30_0 <- tests[tests$score == "MetA_minus_MetB" & tests$comparison == "Dox_30 vs Dox_0", ]
comparison <- data.frame(
  metric = c("Signature genes found on array",
             "MetA - MetB score, 30 vs 0 ng/mL: direction",
             "MetA - MetB score, 30 vs 0 ng/mL: p-value",
             "Highest-MetA group"),
  ours  = c(paste(nrow(used), "of", nrow(signature)),
            ifelse(diff_30_0$difference > 0, "towards MetA", "towards MetB"),
            signif(diff_30_0$p_value, 3),
            GROUP_LEVELS[which.max(tapply(scores$MetA_minus_MetB, scores$group, mean))]),
  paper = c("157", "towards MetA (~60% MetB -> ~70% MetA)", "< 0.001", "Dox_30")
)
print(comparison, row.names = FALSE)
write.csv(comparison, file.path(DIR_RESULTS, "aim2_comparison_with_paper.csv"), row.names = FALSE)

# ---- 6. Figures ----------------------------------------------------------------------------
plot_df <- do.call(rbind, lapply(intersect(c("MetA", "MetB"), subtypes), function(s)
  data.frame(sample = scores$sample, group = scores$group, subtype = s, score = scores[[s]])))
means_df <- aggregate(score ~ group + subtype, plot_df, mean)
p_scores <- ggplot(plot_df, aes(group, score, fill = subtype)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_col(data = means_df, position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             size = 1.8, show.legend = FALSE) +
  scale_fill_manual(values = c(MetA = "#2E75B6", MetB = "#C00000")) +
  scale_x_discrete(labels = GROUP_LABELS) +
  labs(x = NULL, y = "Subtype score (mean z-score of signature genes)", fill = NULL,
       title = "MetA and MetB signature scores (cf. paper Fig. 3C)",
       subtitle = paste0(nrow(used), " of ", nrow(signature),
                         " Thysell et al. (2022) genes; scores are relative to these 16 arrays")) +
  theme_bw()
ggsave(file.path(DIR_FIGURES, "aim2_subtype_scores.png"), p_scores, width = 8, height = 5, dpi = 150)

ht <- Heatmap(z, name = "z-score",
              col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
              row_split = used$subtype, cluster_row_slices = FALSE,
              show_row_names = nrow(z) <= 80, row_names_gp = gpar(fontsize = 6),
              column_split = eset$group, cluster_column_slices = FALSE,
              top_annotation = HeatmapAnnotation(Group = eset$group,
                                                 col = list(Group = GROUP_COLORS)),
              column_title = "MetA-C signature genes (Thysell et al., 2022)")
png(file.path(DIR_FIGURES, "aim2_signature_heatmap.png"), width = 1500, height = 1700, res = 150)
draw(ht)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(DIR_RESULTS, "sessionInfo_04.txt"))
message("Aim 2 finished. Tables in ", DIR_RESULTS, "/, figures in ", DIR_FIGURES, "/")
