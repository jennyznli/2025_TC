# ============================================================
#   This script performs joint processing/imputation for
#   pediatric + adult betas, and creates a joint spreadsheet.
# ============================================================
source("../config.R")

R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

# ========================
# LOAD DATA
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

# ========================
# MAKE JOINT SPREADSHEET
# ========================
cols <- c("Source","Batch", "IDAT", "Sample_ID", "Lymph_Node", "Paired_Primary", 
          "Consensus_Cluster", "Chronological_Age", "Sex",
          "Histology", "Clinical_Invasiveness", "N", "M", "Driver_Group",
          "Epigenetic_Age", "IC_EpiDISH", "Probe_Success_Rate", "Reference")

ss_adt <- read_excel(file.path("ss", ADT_META))
ss_adt$Consensus_Cluster <- NA
ss_adt$Batch <- "ADT"
ss_adt <- ss_adt[, cols]

ss_ped <- read_excel(file.path("ss", PED_META))
ss_ped <- ss_ped[, cols]
dim(ss_ped)
ss_ped$Chronological_Age <- as.numeric(ss_ped$Chronological_Age)

ss <- rbind(ss_ped, ss_adt)

write_xlsx(ss, file.path("ss", "joint_master.xlsx"))

