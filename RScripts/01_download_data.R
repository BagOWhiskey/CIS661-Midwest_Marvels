# =============================================================================
# 01_download_data.R
# Downloads the raw CEL files and sample metadata for GSE272451 from GEO and
# builds data/sample_sheet.csv (one row per array: GSM ID, group, replicate,
# CEL file path).
#
# Download size: ~420 MB. Run once; later runs reuse the files on disk.
# =============================================================================

# Work from the project folder even if RStudio has the parent project open
if (!file.exists(file.path("R", "_config.R")) &&
    file.exists(file.path("mir375_replication", "R", "_config.R"))) setwd("mir375_replication")

source(file.path("R", "_config.R"))
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

options(timeout = 3600)  # large download; R's default 60-second timeout is too short

# ---- 1. Raw CEL files ------------------------------------------------------------
tar_file <- file.path(DIR_RAW, GSE_ID, paste0(GSE_ID, "_RAW.tar"))
if (!file.exists(tar_file)) {
  getGEOSuppFiles(GSE_ID, baseDir = DIR_RAW)   # the series' only file is the RAW.tar
}
stopifnot(file.exists(tar_file))

cel_files <- list.files(DIR_CEL, pattern = "\\.cel(\\.gz)?$", ignore.case = TRUE)
if (length(cel_files) == 0) {
  untar(tar_file, exdir = DIR_CEL)
  cel_files <- list.files(DIR_CEL, pattern = "\\.cel(\\.gz)?$", ignore.case = TRUE)
}
message(length(cel_files), " CEL files in ", DIR_CEL)

# ---- 2. Sample metadata (series matrix) ------------------------------------------------
gse   <- getGEO(GSE_ID, GSEMatrix = TRUE, destdir = DIR_RAW)[[1]]
pdata <- pData(gse)

# Keep the raw GEO metadata for the record (and for checking the parsing below)
write.csv(pdata, file.path(DIR_RAW, paste0(GSE_ID, "_geo_metadata.csv")))

# ---- 3. Build the sample sheet from the sample titles -----------------------------
# Titles look like "C4-2B, WT_1" or "C4-2B, Dox_30_4"
parse_title <- function(title) {
  label     <- trimws(sub("^.*,", "", title))   # "WT_1" / "Dox_30_4"
  replicate <- as.integer(sub("^.*_", "", label))
  group     <- sub("_[0-9]+$", "", label)       # "WT" / "Dox_30"
  data.frame(group = group, replicate = replicate)
}

sample_sheet <- data.frame(
  gsm   = pdata$geo_accession,
  title = pdata$title,
  parse_title(pdata$title),
  stringsAsFactors = FALSE
)

# Match each GSM to its CEL file (GEO names the files "GSMxxxxxxx_<name>.CEL.gz")
cel_gsm <- regmatches(cel_files, regexpr("GSM[0-9]+", cel_files))
sample_sheet$cel_file <- file.path(DIR_CEL, cel_files[match(sample_sheet$gsm, cel_gsm)])

# Short, readable sample IDs for plots: WT_1, Dox0_1, ...
sample_sheet$sample_id <- paste0(sub("_", "", sample_sheet$group), "_", sample_sheet$replicate)
sample_sheet$group     <- factor(sample_sheet$group, levels = GROUP_LEVELS)
sample_sheet           <- sample_sheet[order(sample_sheet$group, sample_sheet$replicate), ]

# ---- 4. Sanity checks: stop early if anything is off -----------------------------------
if (any(is.na(sample_sheet$group))) {
  print(sample_sheet[is.na(sample_sheet$group), c("gsm", "title")])
  stop("Some sample titles did not match GROUP_LEVELS; check parse_title().")
}
if (any(is.na(match(sample_sheet$gsm, cel_gsm)))) {
  stop("Some GSM IDs have no CEL file; check the files in ", DIR_CEL)
}
print(table(sample_sheet$group))
stopifnot(nrow(sample_sheet) == 16, all(table(sample_sheet$group) == 4))

write.csv(sample_sheet, SAMPLE_SHEET_FILE, row.names = FALSE)
message("Sample sheet written to ", SAMPLE_SHEET_FILE)
