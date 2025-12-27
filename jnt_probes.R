# ============================================================
#   This script performs joint processing/imputation for
#   pediatric + adult betas.
# ============================================================

# ========================
# CONFIGURATION
# ========================
BASE_DIR <- "/home/lijz/labprojects/20250721_Jenny/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# LOADING DATA
# ========================
betas_ped <- readRDS(file.path(DATA_DIR, "ped_betas_QCDPB_prc.rds"))
cat("Collapsed EPICv2 beta dimensions:", dim(betas_ped), "\n")

betas_adult <- readRDS(file.path(DATA_DIR, "adt_betas_QCDPB_pr.rds"))
cat("Original HM450 beta dimensions:", dim(betas_adult), "\n")

# ========================
# COMBINING BETAS
# ========================
ped_probes <- rownames(betas_ped)
adult_probes <- rownames(betas_adult)
cat("Pediatric probe count:", length(ped_probes), "\n")
cat("Adult probe count:", length(adult_probes), "\n")

joint_probes <- intersect(ped_probes, adult_probes)
saveRDS(joint_probes, file.path(DATA_DIR, "jnt_QCDPB_pr_probes.rds"))
cat("Overlapping probes:", length(joint_probes), "\n")

joint_betas <- cbind(betas_ped[joint_probes,], betas_adult[joint_probes,])
cat("Joint beta matrix dimensions:", dim(joint_betas), "\n")

saveRDS(joint_betas, file.path(DATA_DIR, "jnt_betas_QCDPB_pr.rds"))

