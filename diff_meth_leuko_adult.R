# ============================================================
#   FIG 2.
#   This script performs differential methylation analysis
#   for 88 primary tumor samples:
#   Invasiveness & Methylation Clusters
# ============================================================
# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- "/home/lijz/labprojects/20250721_Jenny/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# INVASIVENESS
# ========================
cat("=== CLINICAL INVASIVENESS ANALYSIS ===\n")

# Load data
ss <- read_excel(file.path(SS_DIR, "adult_ln_master.xlsx"))
ss_primary <- ss %>% filter(Lymph_Node == "F")
betas <- readRDS(file.path(DATA_DIR, "adt_betas_QCDPB_pr.rds"))
print(dim(ss_primary))

ss_primary$Chronological_Age <- as.numeric(ss_primary$Chronological_Age)
ss_primary$Sex <- as.factor(ss_primary$Sex)
ss_primary$Clinical_Invasiveness <- as.factor(ss_primary$Clinical_Invasiveness)
ss_primary$IC_EpiDISH <- as.numeric(ss_primary$IC_EpiDISH)

betas <- betas[, ss_primary$IDAT]
print(dim(betas))

se <- SummarizedExperiment(betas, colData = ss_primary)

se_ok <- (checkLevels(assay(se), colData(se)$Sex) &
              checkLevels(assay(se), colData(se)$Clinical_Invasiveness))

cat("Probes passing QC:", sum(se_ok), "out of", length(se_ok), "\n")
se <- se[se_ok, ]
cat("Final SE dimensions:", dim(se), "\n")

colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
colData(se)$Clinical_Invasiveness <- relevel(factor(colData(se)$Clinical_Invasiveness), "Low")

cat("Reference levels - Sex:", levels(colData(se)$Sex)[1],
    "Clinical_Invasiveness:", levels(colData(se)$Clinical_Invasiveness)[1], "\n")

smry <- DML(se, ~Clinical_Invasiveness + Sex  + IC_EpiDISH,
            BPPARAM = BiocParallel::MulticoreParam(workers = N_WORKERS))
res_inv <- summaryExtractTest(smry)

cat("Results dimensions:", dim(res_inv), "\n")
cat("First few results:\n")
print(head(res_inv, 3))

saveRDS(res_inv, file = file.path(DATA_DIR, "adt_dm_invasiveness_leuko.rds"))
cat("Clinical_Invasiveness differential methylation completed:", dim(res_inv), "\n")

