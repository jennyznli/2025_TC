# ============================================================
#   This script extracts methylation betas and performs preprocessing.
#   Usage: Rscript extract_idats.R <dataset_name> <processing> <mask>
#   Example: Rscript extract_idats.R adt QCDPB TRUE
# ============================================================
library(sesame)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript extract_idats.R <dataset_name> <processing> <mask>\n",
         "Example: Rscript extract_idats.R adt QCDPB TRUE")
}

# CONFIGS
source("../config.R")

DATASET_NAME <- args[1]
PROCESSING <- args[2]
MASK <- args[3]

MASK_LOGICAL <- as.logical(MASK)
if (is.na(MASK_LOGICAL)) {
    stop("Invalid mask value. Please use TRUE or FALSE.")
}

CONFIG <- BiocParallel::SerialParam()

IDAT_DIR <- file.path(BASE_DIR, "idat", DATASET_NAME)

# LOAD BETAS
betas <- openSesame(IDAT_DIR, func = getBetas, prep = PROCESSING, mask = MASK_LOGICAL, BPPARAM = CONFIG)
cat("Successfully loaded betas.\n")

if (is.null(betas)) {
    stop("Failed to extract beta values")
}
cat("Beta matrix dimensions:", dim(betas), "\n")

# STATS
nas <- sum(is.na(betas))
total_values <- prod(dim(betas))
na_rate <- nas / total_values
dim_vals <- dim(betas)

# CREATE FILENAME BASED ON INPUTS
mask_suffix <- ifelse(MASK_LOGICAL, "_mask", "")
output_filename <- paste0(DATASET_NAME, "_betas_", PROCESSING, mask_suffix, ".rds")
output_path <- file.path(DATA_DIR, output_filename)

saveRDS(betas, output_path)
cat("Saved beta values to:", output_filename, "\n")

# ========================
# FINAL SUMMARY
# ========================
cat("Dataset:", DATASET_NAME, "\n")
cat("Samples:", dim_vals[2], "\n")
cat("Probes:", dim_vals[1], "\n")
cat("NAs:", nas, "\n")
cat("Mask rate:", na_rate, "\n")
