# ============================================================
#   This script performs PCA on all betas and the most variable probes.
#   Usage: Rscript pca.R <betas_filename> <number of top variable probes used>
#   Example: Rscript pca.R ped_betas_QCDPBP_prc.rds 30000
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript pca.R <input_rds_file> <number of top variable probes used>")

# CONFIGS 
source("../config.R")

input <- args[1]
n_probes <- as.numeric(args[2])
cat("Input betas file:", input, "\n")
cat("Input probes selected:", n_probes, "\n")

# READ DATA IN 
data <- readRDS(input)
name <- file_path_sans_ext(basename(input))
dir <- dirname(input)
cat("Initial betas dimensions:", dim(data), "\n")

# PCA on all probes
pca_all <- prcomp(t(data), scale. = TRUE)
output_file <- file.path(dir, paste0(name, "_pca.rds"))
saveRDS(pca_all, output_file)
cat("PCA on all probes dimensions:", dim(pca_all$x), "\n")

# PCA on top most variable probes
top_var <- bSubMostVariable(data, n_probes)
cat("After subsetting betas to most variable dimensions:", dim(top_var), "\n")

pca_var <- prcomp(t(top_var), scale. = TRUE)
cat("PCA on most variable probes dimensions:", dim(pca_var$x), "\n")

output_file_var <- file.path(dir, paste0(name, "_", n_probes, "_pca.rds"))
saveRDS(pca_var, output_file_var) 
cat("PCA on most variable probes dimensions:", dim(pca_var$x), "\n")
