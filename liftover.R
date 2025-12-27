# ============================================================
#   This script lifts over betas.
#   Usage: Rscript liftover.R <betas_file> <platform> 
#   Example: Rscript liftover.R ped_betas_QCDPB_mask_prc_imp.rds HM450
# ============================================================
library(sesame)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript liftover.R <betas_file> <platform>\n",
         "Example: Rscript liftover.R ped_betas_QCDPB_mask_prc_imp.rds HM450")
}

BETAS <- args[1]
PLATFORM <- args[2]
CONFIG <- BiocParallel::SerialParam()

# setup directory
BASE_DIR <- "/home/lijz/labprojects/20250721_Jenny/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

betas <- readRDS(file.path(DATA_DIR, BETAS))
betas_lifted <- mLiftOver(betas, PLATFORM)

base_name <- tools::file_path_sans_ext(BETAS)
output_filename <- paste0(base_name, "_lift", PLATFORM, ".rds")

cat("Saving lifted betas to:", file.path(DATA_DIR, output_filename), "\n")
saveRDS(betas_lifted, file.path(DATA_DIR, output_filename))

cat("Completed.\n")
