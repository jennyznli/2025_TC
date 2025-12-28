# ============================================================
# Performs differential methylation analysis for the adult cohort. 
# ============================================================
source("../config.R")

# ========================
# BH CORRECTION
# ========================
res <- readRDS(here("data", "adt_dm_invasiveness_leuko.rds"))
res$Clinical_InvasivenessHigh_BH <- p.adjust(res$Pval_Clinical_InvasivenessHigh, method = "BH")
saveRDS(res, file.path("data", "adt_dm_invasiveness_leuko_bh.rds"))

# ========================
# ENRICHMENT
# ========================
ss <- read_excel(file.path("ss", ADT_META))
res <- readRDS(file.path("data", "adt_dm_invasiveness_leuko_bh.rds"))

# Filter for significant probes
sig_probes <- res %>%
  dplyr::filter(Clinical_InvasivenessHigh_BH < 0.05, abs(Est_Clinical_InvasivenessHigh) > 0.2) %>%
  mutate(Methylation = ifelse(Est_Clinical_InvasivenessHigh > 0.2, "Hyper", "Hypo"))

hyper <- sig_probes %>%
  filter(Methylation == "Hyper")
hypo <- sig_probes %>%
  filter(Methylation == "Hypo")

# LOW vs HIGH — TFBS HYPO
x <- testEnrichment(hypo$Probe_ID, "TFBSconsensus", platform = "HM450", universe = res$Probe_ID)
pdf(file.path("figures", "adt_dm_inv_hypo_tfbs_enrichment_leuko.pdf"),
    width = 3, height = 4, onefile = FALSE)
plotDotFDR(x, n_min = 30, n_max = 30)
dev.off()
saveRDS(x, file.path("data", "adt_dm_inv_hypo_tfbs_enrichment_leuko.rds"))

# # LOW vs HIGH — TFBS HYPER
# y <- testEnrichment(hyper$Probe_ID, "TFBSconsensus", platform = "EPIC", universe = res$Probe_ID)
# pdf(file.path("figures", "adt_dm_inv_hyper_tfbs_enrichment_leuko.pdf"),
#     width = 3, height = 4, onefile = FALSE)
# plotDotFDR(y, n_min = 30, n_max = 30)
# dev.off()
# saveRDS(y, file.path("data", "adt_dm_inv_hyper_tfbs_enrichment_leuko.rds"))
# 
