# =============================================================================
# 03_limma_de.R
# Aim 1: differential expression between the four groups with limma
# (paper: 7,656 transcripts with adj.P < 0.05 between at least two groups;
#  Fig. 3D heatmap; Fig. 4A volcano plot, 30 vs. 0 ng/mL Dox).
#
# Input:  data/processed/eset_rma_annotated.rds (from 02_qc_normalize.R)
# Output: results/de_*.csv, results/aim1_comparison_with_paper.csv,
#         figures/aim1_*.png, data/processed/limma_fit.rds (used by Aim 3)
# =============================================================================

# Work from the project folder even if RStudio has the parent project open
if (!file.exists(file.path("R", "_config.R")) &&
    file.exists(file.path("mir375_replication", "R", "_config.R"))) setwd("mir375_replication")

source(file.path("R", "_config.R"))
suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})
set.seed(SEED)

eset <- readRDS(ESET_FILE)
eset$group <- factor(eset$group, levels = GROUP_LEVELS)
n_start <- nrow(eset)
symbols_on_array <- unique(fData(eset)$SYMBOL)

# ---- 1. Probe filtering (our choice; the paper does not state its filter) ------------
filter_log <- data.frame(step = "RMA output", probesets = n_start)

if (FILTER_UNANNOTATED) {
  eset <- eset[!is.na(fData(eset)$SYMBOL), ]
  filter_log <- rbind(filter_log, data.frame(step = "with gene symbol", probesets = nrow(eset)))
}
if (FILTER_LOW_EXPR) {
  cutoff <- quantile(exprs(eset), LOW_EXPR_QUANTILE)
  expressed <- rowSums(exprs(eset) > cutoff) >= LOW_EXPR_MIN_ARRAYS
  eset <- eset[expressed, ]
  filter_log <- rbind(filter_log, data.frame(
    step = sprintf("expressed (> %.2f log2 in >= %d arrays)", cutoff, LOW_EXPR_MIN_ARRAYS),
    probesets = nrow(eset)))
}
rownames(filter_log) <- NULL
print(filter_log)

# Warn if the filter removed any of the paper's highlighted genes: that would be
# a consequence of OUR filtering choice, not evidence against the paper
removed_by_filter <- setdiff(intersect(HIGHLIGHT_GENES, symbols_on_array), fData(eset)$SYMBOL)
not_on_array      <- setdiff(HIGHLIGHT_GENES, symbols_on_array)
if (length(removed_by_filter) > 0) {
  warning("Removed by the probe filter: ", paste(removed_by_filter, collapse = ", "),
          ". Consider loosening LOW_EXPR_QUANTILE in R/_config.R.")
}
if (length(not_on_array) > 0) {
  warning("No probeset annotated to: ", paste(not_on_array, collapse = ", "))
}
write.csv(filter_log, file.path(DIR_RESULTS, "aim1_probe_filtering.csv"), row.names = FALSE)

# ---- 2. Linear model: one coefficient per group (no intercept) -------------------------
design <- model.matrix(~ 0 + group, data = pData(eset))
colnames(design) <- levels(eset$group)

# All six pairwise comparisons, so "significant between at least two groups"
# matches the paper's definition
contrast_matrix <- makeContrasts(
  Dox0_vs_WT     = Dox_0  - WT,
  Dox10_vs_WT    = Dox_10 - WT,
  Dox30_vs_WT    = Dox_30 - WT,
  Dox10_vs_Dox0  = Dox_10 - Dox_0,
  Dox30_vs_Dox0  = Dox_30 - Dox_0,    # Fig. 4A
  Dox30_vs_Dox10 = Dox_30 - Dox_10,
  levels = design
)

fit  <- lmFit(eset, design)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))
saveRDS(list(fit = fit, fit2 = fit2, design = design,
             contrasts = contrast_matrix), FIT_FILE)

# ---- 3. Count significant probesets -----------------------------------------------------
# Benjamini-Hochberg within each contrast ("separate"), as limma does by default
dt <- decideTests(fit2, method = "separate", adjust.method = "BH", p.value = ALPHA)
per_contrast <- data.frame(
  contrast = colnames(dt),
  up       = colSums(dt ==  1),
  down     = colSums(dt == -1),
  row.names = NULL
)
per_contrast$total <- per_contrast$up + per_contrast$down
print(per_contrast)
write.csv(per_contrast, file.path(DIR_RESULTS, "aim1_de_counts_per_contrast.csv"), row.names = FALSE)

de_any <- rowSums(dt != 0) > 0

# Alternative definition: one moderated F-test across all four groups
f_table <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "none")
de_f    <- f_table$adj.P.Val < ALPHA

# ---- 4. Result tables ----------------------------------------------------------------------
for (cn in colnames(contrast_matrix)) {
  tt <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH")
  write.csv(tt, file.path(DIR_RESULTS, paste0("de_", cn, ".csv")), row.names = FALSE)
}
write.csv(f_table[de_any, ], file.path(DIR_RESULTS, "de_significant_any_contrast.csv"),
          row.names = FALSE)

# The six genes highlighted in the paper (several probesets per gene are possible)
highlight <- do.call(rbind, lapply(c("Dox10_vs_Dox0", "Dox30_vs_Dox0"), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH")
  tt <- tt[tt$SYMBOL %in% HIGHLIGHT_GENES, c("PROBEID", "SYMBOL", "logFC", "adj.P.Val")]
  cbind(contrast = cn, tt)
}))
highlight$down_and_significant <- highlight$logFC < 0 & highlight$adj.P.Val < ALPHA
highlight$note <- ""
missing_genes <- setdiff(HIGHLIGHT_GENES, highlight$SYMBOL)
if (length(missing_genes) > 0) {
  highlight <- rbind(highlight, data.frame(
    contrast = NA, PROBEID = NA, SYMBOL = missing_genes, logFC = NA, adj.P.Val = NA,
    down_and_significant = NA,
    note = ifelse(missing_genes %in% removed_by_filter, "removed by probe filter",
                  "not annotated on array")))
}
highlight <- highlight[order(highlight$SYMBOL, highlight$contrast), ]
print(highlight, row.names = FALSE)
write.csv(highlight, file.path(DIR_RESULTS, "aim1_highlighted_genes.csv"), row.names = FALSE)

# ---- 5. Comparison with the paper (proposal section E4) ------------------------------------
genes_down_30 <- unique(highlight$SYMBOL[which(highlight$contrast == "Dox30_vs_Dox0" &
                                                 highlight$down_and_significant)])
comparison <- data.frame(
  metric = c("Probesets tested",
             "DE in >= 1 pairwise contrast (adj.P < 0.05)",
             "DE by moderated F-test (adj.P < 0.05)",
             "Highlighted genes down at adj.P < 0.05 (30 vs 0 ng/mL)"),
  ours  = c(nrow(eset), sum(de_any), sum(de_f),
            paste0(length(genes_down_30), " of 6",
                   if (length(missing_genes) > 0)
                     paste0(" (not tested: ", paste(missing_genes, collapse = ", "), ")"))),
  paper = c("not reported", PAPER$n_de_any_contrast, "not reported", "6 of 6")
)
print(comparison, row.names = FALSE)
write.csv(comparison, file.path(DIR_RESULTS, "aim1_comparison_with_paper.csv"), row.names = FALSE)

# ---- 6. Volcano plot, 30 vs. 0 ng/mL Dox (Fig. 4A) -----------------------------------------
vt <- topTable(fit2, coef = "Dox30_vs_Dox0", number = Inf, adjust.method = "BH")
vt$significance <- cut(vt$adj.P.Val, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("adj.P < 0.001", "adj.P < 0.01", "adj.P < 0.05", "n.s."))
vt$label <- ifelse(vt$SYMBOL %in% HIGHLIGHT_GENES, vt$SYMBOL, NA)

p_volcano <- ggplot(vt, aes(logFC, -log10(adj.P.Val), color = significance)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("adj.P < 0.001" = "#1F3FBF", "adj.P < 0.01" = "black",
                                "adj.P < 0.05" = "#A50F15", "n.s." = "grey70")) +
  labs(x = "log2 fold change (30 vs. 0 ng/mL Dox)", y = "-log10(adjusted p-value)",
       color = NULL, title = "miR-375 induction in C4-2B-375 cells (cf. paper Fig. 4A)") +
  theme_bw()
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(aes(label = label), color = "black", fontface = "bold",
                             size = 3.5, na.rm = TRUE, max.overlaps = Inf)
} else {
  p_volcano <- p_volcano + geom_text(aes(label = label), color = "black", na.rm = TRUE,
                                     vjust = -0.6, size = 3.2)
}
ggsave(file.path(DIR_FIGURES, "aim1_volcano_30_vs_0.png"), p_volcano,
       width = 7, height = 6, dpi = 150)

# ---- 7. Heatmap of DE probesets with four gene clusters (Fig. 3D) -----------------------------
# Free memory first: clustering ~9,000 genes is the most memory-hungry step in this
# script, and GVSU sessions are capped at 4 GB
rm(f_table, vt, p_volcano, tt, fit)
invisible(gc())

z <- t(scale(t(exprs(eset)[de_any, , drop = FALSE])))    # row z-scores
row_hc   <- hclust(dist(z), method = "complete")
clusters <- cutree(row_hc, k = 4)
# cutree() numbers clusters by first appearance in the data, but the heatmap labels
# its slices 1-4 from top to bottom. Renumber so the CSV matches the figure.
clusters <- match(clusters, unique(clusters[row_hc$order]))
cluster_genes <- data.frame(PROBEID = rownames(z),
                            SYMBOL  = fData(eset)$SYMBOL[match(rownames(z), featureNames(eset))],
                            cluster = clusters)
write.csv(cluster_genes, file.path(DIR_RESULTS, "aim1_heatmap_gene_clusters.csv"), row.names = FALSE)

col_anno <- HeatmapAnnotation(Group = eset$group,
                              col = list(Group = GROUP_COLORS))
ht <- Heatmap(z, name = "z-score",
              col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
              cluster_rows = row_hc, row_split = 4, show_row_names = FALSE,
              cluster_columns = TRUE, column_split = eset$group,
              cluster_column_slices = FALSE,   # keep groups in WT, 0, 10, 30 order
              top_annotation = col_anno, use_raster = TRUE,
              column_title = paste0(sum(de_any), " DE probesets (adj.P < ", ALPHA,
                                    " in >= 1 contrast)"))
# If drawing fails (e.g. memory), delete the half-drawn file instead of leaving it behind
heatmap_file <- file.path(DIR_FIGURES, "aim1_heatmap_de.png")
png(heatmap_file, width = 1500, height = 1800, res = 150)
heatmap_ok <- tryCatch({ draw(ht); TRUE },
                       error = function(e) { message("Heatmap failed: ", conditionMessage(e)); FALSE })
dev.off()
if (!heatmap_ok) {
  unlink(heatmap_file)
  warning("Heatmap not saved. Restart R (Session -> Restart R) and rerun this script.")
}

# ---- 8. Cluster profiles: what each heatmap cluster does across the four groups ------------
# Mean z-score per group for each cluster, so the heatmap pattern can be described in
# words (e.g. "cluster 3 falls with miR-375 dose") and compared with the paper's Fig. 3D
group_means <- sapply(GROUP_LEVELS, function(g) rowMeans(z[, eset$group == g, drop = FALSE]))
profile <- do.call(rbind, lapply(sort(unique(clusters)), function(k) {
  m <- group_means[clusters == k, , drop = FALSE]
  data.frame(cluster = k, n_probesets = nrow(m), group = GROUP_LEVELS,
             mean_z = colMeans(m), sd_z = apply(m, 2, sd))
}))
profile$group <- factor(profile$group, levels = GROUP_LEVELS)
profile$cluster_label <- paste0("Cluster ", profile$cluster, " (n = ", profile$n_probesets, ")")

# Which cluster each of the paper's six highlighted genes falls in
hl_clusters <- cluster_genes[cluster_genes$SYMBOL %in% HIGHLIGHT_GENES, c("SYMBOL", "cluster")]
print(hl_clusters, row.names = FALSE)

cluster_summary <- reshape(profile[, c("cluster", "n_probesets", "group", "mean_z")],
                           idvar = c("cluster", "n_probesets"), timevar = "group",
                           direction = "wide")
names(cluster_summary) <- sub("^mean_z\\.", "mean_z_", names(cluster_summary))
cluster_summary$highlighted_genes <- vapply(cluster_summary$cluster, function(k)
  paste(sort(unique(hl_clusters$SYMBOL[hl_clusters$cluster == k])), collapse = ", "), character(1))
rownames(cluster_summary) <- NULL
print(cluster_summary, digits = 2)
write.csv(cluster_summary, file.path(DIR_RESULTS, "aim1_cluster_summary.csv"), row.names = FALSE)

p_profile <- ggplot(profile, aes(group, mean_z, group = 1)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_errorbar(aes(ymin = mean_z - sd_z, ymax = mean_z + sd_z), width = 0.15, color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(color = group), size = 3) +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = c(WT = "WT", Dox_0 = "0", Dox_10 = "10", Dox_30 = "30")) +
  facet_wrap(~ cluster_label, nrow = 1) +
  labs(x = "Group (Dox, ng/mL)", y = "Mean row z-score (+/- SD)",
       title = "Expression pattern of each heatmap cluster") +
  theme_bw()
ggsave(file.path(DIR_FIGURES, "aim1_cluster_profiles.png"), p_profile,
       width = 10, height = 3.5, dpi = 150)

writeLines(capture.output(sessionInfo()), file.path(DIR_RESULTS, "sessionInfo_03.txt"))
message("Aim 1 finished. Tables in ", DIR_RESULTS, "/, figures in ", DIR_FIGURES, "/")
