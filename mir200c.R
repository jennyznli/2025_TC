
# ============================================================
#   This script calculates MIR200C methylation using CytoMethIC models.
#   Usage: Rscript mir200c.R <betas_filename> <platform>
#   Example: Rscript mir200c.R ped_betas_QCDPB_prc.rds EPICv2
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript mir200c.R <betas_filename> <platform>\n",
       "Example: Rscript mir200c.R ped_betas_QCDPB_prc.rds EPICv2")
}

BETAS <- args[1]
PLATFORM <- args[2]
CONFIG <- BiocParallel::SerialParam()

BASE_DIR <- "/home/lijz/labprojects/20250721_Jenny/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "load_packages.R"))

base_name <- tools::file_path_sans_ext(BETAS)
file_ext <- tools::file_ext(BETAS)
output_filename <- paste0(base_name, "_mir200c.", file_ext)

cat("Processing betas file:", BETAS, "\n")
cat("Input Platform:", PLATFORM, "\n")
cat("Output filename:", output_filename, "\n")

# ========================
# LOAD DATA
# ========================
betas <- readRDS(file.path(DATA_DIR, BETAS))
cat("Loaded betas matrix with", nrow(betas), "probes and", ncol(betas), "samples\n")

# ========================
# LOAD MODEL
# ========================
cat("Loading MIR200C CytoMethIC model...\n")
cmi <- readRDS(url("https://github.com/zhou-lab/CytoMethIC_models/raw/refs/heads/main/models/MIR200C_EPIC_20240315.rds"))
cat("Model loaded successfully.\n")

# ========================
# LIFTOVER IF NEEDED
# ========================
# The MIR200C model is trained on EPIC platform
# If input is not EPIC, need to liftover
if (PLATFORM != "EPIC") {
  cat("Platform is", PLATFORM, "- performing liftover to EPIC...\n")
  betas <- mLiftOver(betas, "EPIC")
  cat("Liftover completed.\n")
} else {
  cat("Platform is already EPIC - no liftover needed.\n")
}

# ========================
# RUN PREDICTION
# ========================
cat("Predicting MIR200C methylation...\n")
mir_results <- cmi_predict(betas, cmi)
saveRDS(mir_results, file.path(DATA_DIR, output_filename))

# ========================
# FORMAT AND SAVE RESULTS
# ========================
# Extract predictions and format as data frame
# The cmi_predict function returns a list with predictions
# Typically the result contains a 'prediction' or 'score' element

save_cytomethic_results <- function(pred_list, filename, sample_names = NULL) {
    df <- as.data.frame(do.call(rbind, lapply(pred_list, as.data.frame)))
    if (!is.null(sample_names)) {
        rownames(df) <- sample_names
    } else {
        rownames(df) <- names(pred_list)
    }
    filepath <- file.path(DATA_DIR, filename)
    write.csv(df, filepath, row.names = TRUE)
    cat("Saved", nrow(df), "Cytomethic predictions to", filename, "\n")
    return(df)
}

mir_df <- save_cytomethic_results(mir_results, output_filename, colnames(betas))

cat("MIR200C predictions completed for", ncol(betas), "samples\n")
print(head(mir_df))

# Save as RDS
saveRDS(mir_df, file.path(DATA_DIR, output_filename))
cat("Results saved to:", file.path(DATA_DIR, output_filename), "\n")

# Also save as CSV for easy viewing
csv_filename <- paste0(base_name, "_mir200c.csv")
write.csv(mir_df, file.path(DATA_DIR, csv_filename), row.names = FALSE)
cat("Results also saved as CSV:", file.path(DATA_DIR, csv_filename), "\n")

# ========================
# SUMMARY STATISTICS
# ========================
cat("\n=== Summary Statistics ===\n")
cat("Mean MIR200C score:", mean(mir_df$mir200c_score, na.rm = TRUE), "\n")
cat("Median MIR200C score:", median(mir_df$mir200c_score, na.rm = TRUE), "\n")
cat("Range:", range(mir_df$mir200c_score, na.rm = TRUE), "\n")
cat("Number of NA values:", sum(is.na(mir_df$mir200c_score)), "\n")


