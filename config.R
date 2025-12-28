# Configuration file for TC 2025

setwd("/Users/lijz/Documents/thyroid")

source(file.path("R", "load_packages.R"))
source(file.path("R", "functions.R"))
source(file.path("R", "color_keys.R"))

# ========================
# DATA FILES
# ========================
PED_META <- "20231102_thyroid_master.xlsx"
ADT_META <- "adult_ln_master.xlsx"
JNT_META <- "joint_master.xlsx"

PED_BETAS <- "ped_betas_QCDPB_prc.rds"
ADT_BETAS <- "ped_betas_QCDPB_prc.rds"

