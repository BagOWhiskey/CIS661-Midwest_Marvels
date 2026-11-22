# =============================================================================
# 00_setup.R
# Installs the packages this project needs (only the ones that are missing)
# and records their versions. Run once per computer / RStudio Server account.
# =============================================================================

# Work from the project folder even if RStudio has the parent project open
if (!file.exists(file.path("R", "_config.R")) &&
    file.exists(file.path("mir375_replication", "R", "_config.R"))) setwd("mir375_replication")

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c("BiocManager", "ggplot2", "ggrepel", "ff",   # ff: low-memory mode in 02
               "readxl")                                     # readxl: Table S2 in 04
bioc_pkgs <- c(
  "Biobase", "GEOquery",               # GEO download + ExpressionSet
  "oligo", "pd.clariom.d.human",       # read Clariom D CEL files + RMA
  "AnnotationDbi",
  "clariomdhumantranscriptcluster.db", # probeset -> gene symbol
  "limma",                             # differential expression
  "ComplexHeatmap", "circlize",        # heatmaps
  "org.Hs.eg.db",                      # gene aliases (04) and GO annotation (05)
  "clusterProfiler"                    # GO enrichment (05); first install takes ~10 min
)

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

# ---- Record versions (needed for the reproducibility requirement) -------------
source(file.path("R", "_config.R"))
all_pkgs <- c(cran_pkgs, bioc_pkgs)
versions <- data.frame(
  package = all_pkgs,
  version = vapply(all_pkgs, function(p) as.character(packageVersion(p)), character(1)),
  row.names = NULL
)
versions <- rbind(
  data.frame(package = c("R", "Bioconductor"),
             version = c(paste(R.version$major, R.version$minor, sep = "."),
                         as.character(BiocManager::version()))),
  versions
)
print(versions)
write.csv(versions, file.path(DIR_RESULTS, "package_versions.csv"), row.names = FALSE)

# The paper reports R 4.2.2 and limma 3.64.1 for the differential expression.
# Those two cannot be combined (limma 3.64 needs a newer R), so note whatever
# versions you actually get here in the lab report.
message("limma version installed: ", packageVersion("limma"),
        " (paper used 3.64.1 for the cell-line DE analysis)")
