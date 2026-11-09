# ================================================================
# One-time installation of the R packages used by both apps.
# Run this once in the RStudio Console:  source("setup/install_packages.R")
# ================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(c("shiny", "DT", "ggplot2", "zip", "rJava", "png"))

BiocManager::install(c("EBImage", "RBioFormats"), ask = FALSE, update = FALSE)

cat("\nDone. Restart R (Session > Restart R) before starting an app.\n")
