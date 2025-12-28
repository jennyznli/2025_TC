# ============================================================
# Performs differential methylation for clusters on 
# the pediatric cohort (LI as reference). 
# ============================================================
source("../config.R")

N_WORKERS <- 20

R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

# ========================
# LOAD DATA
# ========================
ss <- read_excel(file.path(SS_DIR, PED_META))
ss_primary <- ss %>% filter(Lymph_Node == "F", Batch == "REF")
betas <- readRDS(file.path(DATA_DIR, PED_BETAS))

ss_primary$Chronological_Age <- as.numeric(ss_primary$Chronological_Age)
ss_primary$Sex <- as.factor(ss_primary$Sex)
ss_primary$Consensus_Cluster <- as.factor(ss_primary$Consensus_Cluster)
ss_primary$IC_EpiDISH <- as.numeric(ss_primary$IC_EpiDISH)

betas <- betas[, ss_primary$IDAT]
print(dim(betas))

se <- SummarizedExperiment(betas, colData = ss_primary)

se_ok <- (checkLevels(assay(se), colData(se)$Sex) &
              checkLevels(assay(se), colData(se)$Consensus_Cluster))

cat("Probes passing QC:", sum(se_ok), "out of", length(se_ok), "\n")
se <- se[se_ok, ]
cat("Final SE dimensions:", dim(se), "\n")

colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
colData(se)$Consensus_Cluster <- relevel(factor(colData(se)$Consensus_Cluster), "LI")

smry <- DML(se, ~Consensus_Cluster + Sex,
            BPPARAM = BiocParallel::MulticoreParam(workers = N_WORKERS))
res <- summaryExtractTest(smry)

cat("Results dimensions:", dim(res), "\n")
cat("First few results:\n")
print(head(res, 3))

saveRDS(res, file = file.path(DATA_DIR, "ped_dm_cluster_li.rds"))

