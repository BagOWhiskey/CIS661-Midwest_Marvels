# =============================================================================
# _config.R
# Shared paths and analysis settings for the miR-375 replication project.
# Every numbered script sources this file first, so change settings HERE only.
#
# Replicating: Järemo et al. (2026), Scientific Reports 16:26685,
#              https://doi.org/10.1038/s41598-026-67032-9
# Dataset:     GEO GSE272451 (C4-2B WT vs. C4-2B-375 at 0/10/30 ng/mL Dox)
# =============================================================================

# ---- Dataset -----------------------------------------------------------------
GSE_ID <- "GSE272451"

# Group order used everywhere (labels come from the GEO sample titles,
# e.g. "C4-2B, WT_1", "C4-2B, Dox_0_1", "C4-2B, Dox_10_1", "C4-2B, Dox_30_1").
GROUP_LEVELS <- c("WT", "Dox_0", "Dox_10", "Dox_30")
GROUP_LABELS <- c(WT = "Wild type", Dox_0 = "0 ng/mL Dox",
                  Dox_10 = "10 ng/mL Dox", Dox_30 = "30 ng/mL Dox")
GROUP_COLORS <- c(WT = "#7F7F7F", Dox_0 = "#1F3FBF",
                  Dox_10 = "#7A1F4F", Dox_30 = "#D62728")   # close to the paper's Fig. 3

# ---- Folder layout (relative to the project root) ------------------------------
DIR_RAW       <- file.path("data", "raw")         # GEO downloads (not committed)
DIR_CEL       <- file.path("data", "raw", "cel")  # extracted CEL files
DIR_PROCESSED <- file.path("data", "processed")   # normalized data (.rds)
DIR_RESULTS   <- "results"                        # tables (.csv) and logs (.txt)
DIR_FIGURES   <- "figures"                        # plots (.png)
for (d in c(DIR_RAW, DIR_CEL, DIR_PROCESSED, DIR_RESULTS, DIR_FIGURES)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

SAMPLE_SHEET_FILE <- file.path("data", "sample_sheet.csv")
ESET_FILE         <- file.path(DIR_PROCESSED, "eset_rma_annotated.rds")
FIT_FILE          <- file.path(DIR_PROCESSED, "limma_fit.rds")

# ---- Memory ----------------------------------------------------------------------
# GVSU RStudio sessions are capped at 4 GB, and reading 16 Clariom D arrays comes
# close to that. 02_qc_normalize.R ran successfully with this set to FALSE (in a fresh
# session with nothing else loaded). If it ever stops with "Memory Limit Exceeded",
# set this to TRUE: the ff package then keeps the probe-level data on disk and oligo
# processes 4 arrays at a time.
LOW_MEMORY_MODE <- FALSE
FF_DIR          <- file.path(DIR_PROCESSED, "ff_tmp")   # temporary, deleted after RMA

# ---- Analysis settings ----------------------------------------------------------
SEED      <- 661      # set.seed() wherever randomness is involved
ALPHA     <- 0.05     # Benjamini-Hochberg adjusted p-value cutoff (as in the paper)

# Probe filtering. The paper does not say how it filtered the cell-line arrays,
# so these are OUR choices: report them as a deviation and try other values as
# a sensitivity analysis (proposal, Section F).
FILTER_UNANNOTATED  <- TRUE  # drop probesets with no gene symbol (mostly controls)
FILTER_LOW_EXPR     <- TRUE  # drop probesets that are never clearly expressed
LOW_EXPR_QUANTILE   <- 0.20  # "expressed" = above this quantile of all log2 values...
LOW_EXPR_MIN_ARRAYS <- 4     # ...in at least this many arrays (= one full group)

# Genes highlighted in the paper (Fig. 4A, Table S5)
HIGHLIGHT_GENES <- c("ACSL3", "YAP1", "RBM24", "RMDN3", "QKI", "CHSY1")

# Numbers reported in the paper, used for the comparison table
PAPER <- list(
  n_de_any_contrast   = 7656,  # "adj.P < 0.05 between at least two groups"
  n_pred_targets      = 1058,  # predicted miR-375 targets among DE transcripts (Aim 3)
  n_stringent_targets = 120
)

# ---- Aim 2: MetA-C signature (04_subtype_metab.R) --------------------------------
# Table S2 of Thysell et al. (2022), Mol Oncol 16:846-859, "Gene transcripts used for
# classification of the metastasis subtypes A-C" (157 transcripts). Download it from
# the article's Supporting Information (https://doi.org/10.1002/1878-0261.13158) and
# save it as this file (.xlsx as downloaded, or saved as .csv).
SIGNATURE_FILE <- file.path("data", "Thysell2022_TableS2.xlsx")
# Leave NULL to detect the columns automatically; otherwise give the exact column
# names in Table S2 that hold the gene symbol and the subtype (MetA/MetB/MetC)
SIGNATURE_GENE_COL    <- NULL
SIGNATURE_SUBTYPE_COL <- NULL

# ---- Aim 3: miR-375 targets (05_targets_go.R) ------------------------------------
# TargetScanHuman 8.0 summary file (downloaded automatically if the server can reach
# targetscan.org; otherwise download it from https://www.targetscan.org/vert_80/ ->
# "Data download" and save it under this name)
TARGETSCAN_URL  <- "https://www.targetscan.org/vert_80/vert_80_data_download/Summary_Counts.default_predictions.txt.zip"
TARGETSCAN_FILE <- file.path(DIR_RAW, "Summary_Counts.default_predictions.txt.zip")
MIR_NAME        <- "hsa-miR-375-3p"
