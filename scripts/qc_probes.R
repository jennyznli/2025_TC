# ============================================================
#   This script performs probe preprocessing for betas.
#   Usage: Rscript qc_probes.R <betas_filename> <platform>
#   Example: Rscript qc_probes.R ped_betas_QCDPB.rds EPICv2
# ============================================================
library(sesame)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript qc_betas.R <betas_filename> <platform>\n",
         "Example: Rscript qc_betas.R ped_betas_QCDPB.rds EPICv2")
}

# CONFIGS
source("../config.R")

BETAS_FILENAME <- args[1]
PLATFORM <- args[2]

CONFIG <- BiocParallel::SerialParam()

R_DIR <- file.path(BASE_DIR, "R")
DATA_DIR <- file.path(BASE_DIR, "data")

# FILENAME
base_name <- tools::file_path_sans_ext(BETAS_FILENAME)
file_ext <- tools::file_ext(BETAS_FILENAME)
output_filename <- paste0(base_name, "_pr.", file_ext)

cat("Input filename:", BETAS_FILENAME, "\n")
cat("Output filename:", output_filename, "\n")

# LOADING DATA
mft <- sesameAnno_readManifestTSV(paste0(PLATFORM, ".hg38.manifest"))

# CONSTRUCT INPUT
input_file <- file.path(DATA_DIR, BETAS_FILENAME)
betas <- readRDS(input_file)
cat("Loaded betas from:", input_file, "\n")
cat("Initial beta dimensions:", dim(betas), "\n")

# PREPROCESSING
# get sex chromosome probes
cat("1. Removing sex chromosome probes...\n")
sex <- mft %>%
    filter(CpG_chrm %in% c("chrX", "chrY"))
sex_probes <- sex$Probe_ID
cat("Found", length(sex_probes), "sex chromosome probes.\n")

betas <- betas[!rownames(betas) %in% sex_probes, ]
cat("After removing sex chromosome probes, beta dimensions:", dim(betas), "\n")

# probe type analysis
cat("\n2. Analyzing probe types...\n")
probes <- rownames(betas)
probe_types <- table(substr(probes, 1, 2))
cat("Probe types before filtering:\n")
print(probe_types)

# Filter to only CG probes
cat("\n3. Filtering to CG probes only...\n")
cg_probes <- probes[substr(probes, 1, 2) == "cg"]
betas <- betas[cg_probes, ]
cat("After keeping only CG probes, beta dimensions:", dim(betas), "\n")
cat("Kept", length(cg_probes), "CG probes out of", length(probes), "total probes\n")

# EPICv2 PROCESSING
if (PLATFORM == "EPICv2") {
    cat("\n4. EPICv2 detected, collapsing to probe prefixes ...\n")
    betas_collapsed <- betasCollapseToPfx(betas)
    cat("After collapsing, beta dimensions:", dim(betas_collapsed), "\n")
    
    # Create collapsed filename with betasc suffix
    collapsed_base_name <- tools::file_path_sans_ext(output_filename)
    collapsed_filename <- paste0(collapsed_base_name, "c.", file_ext)
    collapsed_output_file <- file.path(DATA_DIR, collapsed_filename)
    
    saveRDS(betas_collapsed, collapsed_output_file)
    cat("Saved collapsed betas to:", collapsed_output_file, "\n")
    
    final_betas <- betas_collapsed
    final_output_file <- collapsed_output_file
    final_filename <- collapsed_filename
} else {
    final_betas <- betas
    final_output_file <- file.path(DATA_DIR, output_filename)
    final_filename <- output_filename
}

cat("\n", ifelse(PLATFORM == "EPICv2", "5", "4"), ". Saving processed data...\n")
output_file <- file.path(DATA_DIR, output_filename)
saveRDS(betas, output_file)
cat("Saved processed betas to:", output_file, "\n")

# SUMMARY
cat("Input filename:", BETAS_FILENAME, "\n")
cat("Platform:", PLATFORM, "\n")
cat("Main output file:", output_filename, "\n")

if (PLATFORM == "EPICv2") {
    cat("Collapsed output file:", final_filename, "\n")
    cat("Final beta dimensions (collapsed):", dim(final_betas), "\n")
} else {
    cat("Final beta dimensions:", dim(betas), "\n")
}

cat("Sex chromosome probes removed:", length(sex_probes), "\n")
cat("Non-CG probes removed:", length(probes) - length(cg_probes), "\n")

if (PLATFORM == "EPICv2") {
    cat("EPICv2 collapse applied: YES\n")
} else {
    cat("EPICv2 collapse applied: NO\n")
}

cat("\nProcessing complete\n")
