# =============================================================================
# 05_targets_go.R
# Aim 3: (A) which biological processes characterise each heatmap gene cluster
#            (paper Fig. 3D: reduced cell cycle in clusters 1-2, increased translation
#            in cluster 3, increased adhesion in cluster 4), and
#        (B) are predicted miR-375 targets among the down-regulated genes, as expected
#            for a miRNA? (paper: 1,058 predicted and 120 stringently predicted targets
#            among the 7,656 DE transcripts; Table S5, Fig. 4A)
#
# Input:  data/processed/limma_fit.rds and results/aim1_heatmap_gene_clusters.csv
#         (both from 03_limma_de.R; rerun 03 first so cluster numbers match the heatmap)
# Output: results/aim3_*.csv, figures/aim3_go_clusters.png, figures/aim3_targets_ecdf.png
#
# Deviations from the paper (report these):
#  - Pathway analysis uses Gene Ontology (biological process) with clusterProfiler,
#    not the commercial MetaCore "process networks", so term names will differ.
#  - The paper does not define "stringently predicted". Here: predicted = any
#    miR-375-3p site in TargetScanHuman 8.0 default predictions; stringent = at least
#    one CONSERVED site.
# =============================================================================

# Work from the project folder even if RStudio has the parent project open
if (!file.exists(file.path("R", "_config.R")) &&
    file.exists(file.path("mir375_replication", "R", "_config.R"))) setwd("mir375_replication")

source(file.path("R", "_config.R"))
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  stop("clusterProfiler is not installed. Run R/00_setup.R first.")
}
suppressPackageStartupMessages({
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})
set.seed(SEED)

fit_obj  <- readRDS(FIT_FILE)
fit2     <- fit_obj$fit2
tested   <- fit2$genes                       # probesets tested in Aim 1 (after filtering)
clusters <- read.csv(file.path(DIR_RESULTS, "aim1_heatmap_gene_clusters.csv"),
                     stringsAsFactors = FALSE)
clusters$ENTREZID <- tested$ENTREZID[match(clusters$PROBEID, tested$PROBEID)]

# Direction of each cluster (from 03's cluster summary, if available)
summary_file <- file.path(DIR_RESULTS, "aim1_cluster_summary.csv")
cluster_info <- data.frame(cluster = sort(unique(clusters$cluster)))
cluster_info$n_probesets <- as.vector(table(clusters$cluster)[as.character(cluster_info$cluster)])
if (file.exists(summary_file)) {
  cs <- read.csv(summary_file)
  change <- cs$mean_z_Dox_30 - cs$mean_z_Dox_0
  cluster_info$pattern <- ifelse(change[match(cluster_info$cluster, cs$cluster)] < 0,
                                 "down with miR-375 (30 vs 0 ng/mL)",
                                 "up with miR-375 (30 vs 0 ng/mL)")
} else {
  cluster_info$pattern <- "rerun 03_limma_de.R for cluster patterns"
}

# =============================================================================
# Part A: GO biological-process enrichment per heatmap cluster
# =============================================================================
universe <- unique(na.omit(as.character(tested$ENTREZID)))
message("GO universe: ", length(universe), " genes tested in Aim 1")

go_tables <- lapply(cluster_info$cluster, function(k) {
  genes <- unique(na.omit(as.character(clusters$ENTREZID[clusters$cluster == k])))
  message("Cluster ", k, ": ", length(genes), " genes -> enrichGO (this can take a minute)")
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = ALPHA, qvalueCutoff = 0.2,
                  minGSSize = 15, maxGSSize = 500, readable = TRUE)
  res <- if (is.null(ego)) data.frame() else as.data.frame(ego)
  if (nrow(res) == 0) return(NULL)
  res <- res[res$p.adjust < ALPHA, c("ID", "Description", "GeneRatio", "BgRatio",
                                     "pvalue", "p.adjust", "Count", "geneID")]
  ratio <- function(x) vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), 1)
  res$fold_enrichment <- ratio(res$GeneRatio) / ratio(res$BgRatio)
  cbind(cluster = k, res[order(res$p.adjust), ])
})
go_all <- do.call(rbind, go_tables)
if (is.null(go_all)) go_all <- data.frame()
write.csv(go_all, file.path(DIR_RESULTS, "aim3_go_enrichment_by_cluster.csv"), row.names = FALSE)

# Top 5 terms per cluster: the analogue of the "top five process networks" in Fig. 3D
top_go <- if (nrow(go_all) > 0) do.call(rbind, lapply(split(go_all, go_all$cluster), head, 5)) else go_all
cluster_info$n_go_terms <- vapply(cluster_info$cluster, function(k) sum(go_all$cluster == k), 1)
cluster_info$top_terms  <- vapply(cluster_info$cluster, function(k)
  paste(head(go_all$Description[go_all$cluster == k], 5), collapse = "; "), "")
print(cluster_info[, c("cluster", "n_probesets", "pattern", "n_go_terms")], row.names = FALSE)
write.csv(cluster_info, file.path(DIR_RESULTS, "aim3_cluster_top_go_terms.csv"), row.names = FALSE)

if (nrow(top_go) > 0) {
  top_go$cluster_label <- factor(paste0("Cluster ", top_go$cluster),
                                 levels = paste0("Cluster ", cluster_info$cluster))
  # p-values can underflow to 0 for very strong enrichments; cap so the colour scale works
  top_go$neg_log10_padj <- -log10(pmax(top_go$p.adjust, 1e-300))
  top_go$term <- ifelse(nchar(top_go$Description) > 55,
                        paste0(substr(top_go$Description, 1, 52), "..."), top_go$Description)
  top_go$term <- factor(top_go$term, levels = rev(unique(top_go$term)))
  p_go <- ggplot(top_go, aes(cluster_label, term, size = Count, color = neg_log10_padj)) +
    geom_point() +
    scale_x_discrete(drop = FALSE) +
    scale_color_gradient(low = "#9ECAE1", high = "#08306B") +
    labs(x = NULL, y = NULL, size = "Genes", color = "-log10(adj. p)",
         title = "Top GO biological processes per heatmap cluster (cf. paper Fig. 3D)") +
    theme_bw() + theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(DIR_FIGURES, "aim3_go_clusters.png"), p_go,
         width = 9, height = max(4, 0.28 * nrow(top_go) + 1.5), dpi = 150)
} else {
  warning("No GO term reached adj.P < ", ALPHA, " in any cluster.")
}

# =============================================================================
# Part B: predicted miR-375 targets (TargetScanHuman 8.0)
# =============================================================================
options(timeout = 1800)
if (!file.exists(TARGETSCAN_FILE)) {
  ok <- tryCatch({ download.file(TARGETSCAN_URL, TARGETSCAN_FILE, mode = "wb"); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) unlink(TARGETSCAN_FILE)
}

if (!file.exists(TARGETSCAN_FILE)) {
  warning("TargetScan file not available, so Part B was skipped. Download ",
          basename(TARGETSCAN_FILE), " from https://www.targetscan.org/vert_80/ (Data download), ",
          "save it as ", TARGETSCAN_FILE, " and rerun.")
} else {
  # The file is large: read it in chunks and keep only the miR-375 rows
  inner   <- unzip(TARGETSCAN_FILE, list = TRUE)$Name[1]
  tmp_dir <- tempfile("targetscan_")
  unzip(TARGETSCAN_FILE, files = inner, exdir = tmp_dir)
  con <- file(file.path(tmp_dir, inner), open = "r")
  header <- readLines(con, n = 1)
  keep <- character(0)
  repeat {
    chunk <- readLines(con, n = 200000)
    if (length(chunk) == 0) break
    keep <- c(keep, chunk[grepl("miR-375", chunk, fixed = TRUE)])
  }
  close(con)
  unlink(tmp_dir, recursive = TRUE)
  ts <- read.delim(text = c(header, keep), check.names = FALSE, stringsAsFactors = FALSE)

  col_like <- function(pattern) {
    hit <- grep(pattern, names(ts), ignore.case = TRUE, value = TRUE)
    if (length(hit) == 0) stop("TargetScan column matching '", pattern, "' not found. Columns: ",
                               paste(names(ts), collapse = " | "))
    hit[1]
  }
  c_sym  <- col_like("^gene symbol$")
  c_cons <- col_like("^total num conserved sites$")
  c_non  <- col_like("^total num nonconserved sites$")
  c_rep  <- col_like("^representative mirna$")
  c_ctx  <- col_like("cumulative weighted context")
  if (any(grepl("^species id$", names(ts), ignore.case = TRUE))) {
    ts <- ts[ts[[col_like("^species id$")]] %in% c(9606, "9606"), ]
  }
  ts <- ts[ts[[c_rep]] == MIR_NAME, ]
  if (nrow(ts) == 0) stop("No TargetScan rows for ", MIR_NAME, "; check MIR_NAME in _config.R")

  # One row per gene (TargetScan lists several transcripts per gene)
  ts_gene <- do.call(rbind, lapply(split(ts, ts[[c_sym]]), function(d) data.frame(
    SYMBOL          = d[[c_sym]][1],
    conserved_sites = max(as.numeric(d[[c_cons]]), na.rm = TRUE),
    all_sites       = max(as.numeric(d[[c_cons]]) + as.numeric(d[[c_non]]), na.rm = TRUE),
    context_score   = suppressWarnings(min(as.numeric(d[[c_ctx]]), na.rm = TRUE)))))
  ts_gene$predicted <- ts_gene$all_sites > 0
  ts_gene$stringent <- ts_gene$conserved_sites > 0
  message(sum(ts_gene$predicted), " predicted and ", sum(ts_gene$stringent),
          " stringently predicted (conserved-site) ", MIR_NAME, " target genes in TargetScan 8.0")

  # Gene-level results for 30 vs 0 ng/mL (best probeset per gene = smallest p-value)
  tt <- topTable(fit2, coef = "Dox30_vs_Dox0", number = Inf, adjust.method = "BH")
  tt <- tt[!is.na(tt$SYMBOL), ]
  tt <- tt[order(tt$P.Value), ]
  tt <- tt[!duplicated(tt$SYMBOL), ]
  tt$predicted <- tt$SYMBOL %in% ts_gene$SYMBOL[ts_gene$predicted]
  tt$stringent <- tt$SYMBOL %in% ts_gene$SYMBOL[ts_gene$stringent]
  tt$class <- ifelse(tt$stringent, "Conserved-site target",
                     ifelse(tt$predicted, "Non-conserved-site target", "Not a predicted target"))

  # DE in at least one contrast (gene level), as in Aim 1
  dt <- decideTests(fit2, method = "separate", adjust.method = "BH", p.value = ALPHA)
  de_sym <- unique(fit2$genes$SYMBOL[rowSums(dt != 0) > 0])
  de_sym <- de_sym[!is.na(de_sym)]

  # Are predicted targets over-represented among genes DOWN at 30 vs 0? (Fisher's exact test)
  sig  <- tt$adj.P.Val < ALPHA
  fisher_row <- function(target_col, direction) {
    moved <- sig & (if (direction == "down") tt$logFC < 0 else tt$logFC > 0)
    tab <- table(factor(tt[[target_col]], c(TRUE, FALSE)), factor(moved, c(TRUE, FALSE)))
    ft <- fisher.test(tab)
    data.frame(targets = ifelse(target_col == "stringent", "conserved-site", "any site"),
               direction = direction,
               targets_in_set = tab["TRUE", "TRUE"], set_size = sum(moved),
               pct_targets_in_set = round(100 * tab["TRUE", "TRUE"] / sum(moved), 1),
               pct_targets_background = round(100 * mean(tt[[target_col]]), 1),
               odds_ratio = unname(ft$estimate), p_value = ft$p.value)
  }
  enrichment <- rbind(fisher_row("predicted", "down"), fisher_row("predicted", "up"),
                      fisher_row("stringent", "down"), fisher_row("stringent", "up"))
  print(enrichment, digits = 3, row.names = FALSE)
  write.csv(enrichment, file.path(DIR_RESULTS, "aim3_target_enrichment_30_vs_0.csv"), row.names = FALSE)

  # Comparison with the paper
  hl <- tt[tt$SYMBOL %in% HIGHLIGHT_GENES, c("SYMBOL", "logFC", "adj.P.Val", "predicted", "stringent")]
  comparison <- data.frame(
    metric = c("Predicted targets among DE genes",
               "Stringently predicted targets among DE genes",
               "Highlighted genes that are conserved-site targets",
               "Targets enriched among genes down at 30 vs 0 (Fisher p)"),
    ours  = c(sum(de_sym %in% ts_gene$SYMBOL[ts_gene$predicted]),
              sum(de_sym %in% ts_gene$SYMBOL[ts_gene$stringent]),
              paste0(sum(hl$stringent), " of ", nrow(hl)),
              signif(enrichment$p_value[enrichment$targets == "any site" &
                                          enrichment$direction == "down"], 3)),
    paper = c(PAPER$n_pred_targets, PAPER$n_stringent_targets, "6 of 6", "not reported")
  )
  print(comparison, row.names = FALSE)
  write.csv(comparison, file.path(DIR_RESULTS, "aim3_comparison_with_paper.csv"), row.names = FALSE)

  top_targets <- tt[tt$stringent & tt$logFC < 0 & sig,
                    c("SYMBOL", "GENENAME", "logFC", "adj.P.Val")]
  top_targets <- top_targets[order(top_targets$logFC), ]
  write.csv(top_targets, file.path(DIR_RESULTS, "aim3_downregulated_conserved_targets.csv"),
            row.names = FALSE)
  message(nrow(top_targets), " conserved-site targets are significantly down at 30 vs 0 ng/mL")

  # Classic miRNA check: targets should shift towards negative fold changes
  tt$class <- factor(tt$class, levels = c("Conserved-site target", "Non-conserved-site target",
                                          "Not a predicted target"))
  n_lab <- table(tt$class)
  levels(tt$class) <- paste0(levels(tt$class), " (n = ", n_lab, ")")
  p_ecdf <- ggplot(tt, aes(logFC, color = class)) +
    stat_ecdf(linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
    coord_cartesian(xlim = c(-1.5, 1.5)) +
    scale_color_manual(values = c("#C00000", "#E69F00", "grey50")) +
    labs(x = "log2 fold change (30 vs. 0 ng/mL Dox)", y = "Cumulative fraction of genes",
         color = NULL, title = paste(MIR_NAME, "predicted targets after miR-375 induction")) +
    theme_bw() + theme(legend.position = "bottom", legend.direction = "vertical")
  ggsave(file.path(DIR_FIGURES, "aim3_targets_ecdf.png"), p_ecdf, width = 6.5, height = 5.5, dpi = 150)
}

writeLines(capture.output(sessionInfo()), file.path(DIR_RESULTS, "sessionInfo_05.txt"))
message("Aim 3 finished. Tables in ", DIR_RESULTS, "/, figures in ", DIR_FIGURES, "/")
