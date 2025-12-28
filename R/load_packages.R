# ============================================================
# Install and load all required CRAN and Bioconductor packages
# ============================================================
cat("Loading required packages...\n")

set.seed(123)
DATE <- format(Sys.Date(), "%Y%m%d")

packages <- list(
    cran = c(
        "dplyr", "ggplot2", "readr", "tidyr", "plotly", "readxl", "here", "stringr",
        "Rtsne", "pheatmap", "RColorBrewer", "ggpubr", "ggrepel", "clusterProfiler",
        "msigdbr", "fgsea", "enrichplot", "pvclust", "tidyverse", "randomForest",
        "RSNNS", "e1071", "caret", "ISLR", "pROC", "shapviz", "future",
        "future.apply", "logger", "treeshap", "viridis", "pals", "sva", "uwot"
    ),
    bioc = c(
        "BiocParallel", "sesame", "DESeq2", "limma", "AnnotationHub", "CytoMethIC",
        "SummarizedExperiment", "knowYourCG", "EpiDISH", "methylclockData", "methylclock", "ExperimentHub"
    )
)

# Function to check and load packages
load_packages <- function(pkg_list) {
    missing_cran <- pkg_list$cran[!sapply(pkg_list$cran, requireNamespace, quietly = TRUE)]
    missing_bioc <- pkg_list$bioc[!sapply(pkg_list$bioc, requireNamespace, quietly = TRUE)]

    # Install missing CRAN packages
    if (length(missing_cran) > 0) {
        install.packages(missing_cran, quiet = TRUE)
        cat("Installed missing CRAN packages: ", paste(missing_cran, collapse = ", "), "\n")
    }

    # Install missing Bioconductor packages
    if (length(missing_bioc) > 0) {
        BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, quiet = TRUE)
        cat("Installed missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "), "\n")
    }

    # Load all packages and print any that failed to load
    all_pkgs <- c(pkg_list$cran, pkg_list$bioc)
    failed_pkgs <- c()

    invisible(sapply(all_pkgs, function(x) {
        if (!requireNamespace(x, quietly = TRUE)) {
            failed_pkgs <<- c(failed_pkgs, x)
        } else {
            suppressPackageStartupMessages(library(x, character.only = TRUE))
        }
    }))

    if (length(failed_pkgs) > 0) {
        cat("Failed to load the following packages: ", paste(failed_pkgs, collapse = ", "), "\n")
    } else {
        message("All packages loaded successfully!")
    }
}

load_packages(packages)
