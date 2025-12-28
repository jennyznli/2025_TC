# ============================================================
#   This script gets QC stats.
#   Usage: Rscript stats.R <dataset_name>
# ============================================================
library(sesame)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript stats.R <dataset_name>\nExample: Rscript stats.R adt")
}

DATASET <- args[1]

# CONFIGS
source("../config.R")
DATASET_NAME <- args[1]

CONFIG <- BiocParallel::SerialParam()

DATA_DIR <- file.path(BASE_DIR, "data")
IDAT_DIR <- file.path(BASE_DIR, "idat", DATASET_NAME)

if (!dir.exists(IDAT_DIR)) {
  stop("IDAT directory not found: ", IDAT_DIR)
}

# SESAME QC & STATS
sdf <- openSesame(IDAT_DIR, func = NULL, BPPARAM = CONFIG)
cat("Created SigDF objects.")

stats <- lapply(sdf, sesameQC_calcStats)
cat("Calculated statistics.")

saveRDS(stats, file.path(DATA_DIR, paste0(DATASET, "_stats.rds")))
cat("Saved statistics.\n")

df_list <- lapply(stats, function(item) {
  as.data.frame(sesameQC_getStats(item))
})
combined_df <- do.call(rbind, df_list)

write.csv(combined_df, file.path(DATA_DIR, paste0(DATASET, "_stats.csv")), row.names = TRUE)
cat("Saved statistics CSV.\n")
