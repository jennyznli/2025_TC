# ============================================================
# Cell type deconvolution for pediatric and adult cohorts.
# ============================================================
source("config.R")

# ========================
# LOAD DATA
# ========================
ss_ped <- read_excel(file.path("ss", PED_META))
ss_adt <- read_excel(file.path("ss", ADT_META))

data(centEpiFibIC.m) # solid tissue - epithelial, fibroblast, immune cell types
data(centBloodSub.m) # immune sells - 7 types

betas_ped <- readRDS(file.path("data", "ped_betas_QCDPB_prc.rds"))
betas_adt <- readRDS(file.path("data", "adt_betas_QCDPB_pr.rds"))

### BROAD CELL TYPE DECONVOLUTION ### 
cef_ped <- epidish(beta.m = betas_ped, ref.m = centEpiFibIC.m, method = "RPC")
cef_adt <- epidish(beta.m = betas_adt, ref.m = centEpiFibIC.m, method = "RPC")

# Epi = Epithelial cells - The main parenchymal/functional cells
# Fib = Fibroblasts - Stromal/connective tissue cells
# IC = Immune Cells - All immune/inflammatory cells grouped together

### IMMUNE CELL TYPE DECONVOLUTION ### 
immune_ped <- hepidish(beta.m = betas_ped, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
immune_adt <- hepidish(beta.m = betas_adt, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')

cef_ped_df <- as.data.frame(cef_ped$estF)
cef_ped_df$IDAT <- rownames(cef_ped_df)
cef_adt_df <- as.data.frame(cef_adt$estF)
cef_adt_df$IDAT <- rownames(cef_adt_df)

immune_ped <- as.data.frame(immune_ped)
immune_adt <- as.data.frame(immune_adt)
immune_ped$IDAT <- rownames(immune_ped)
immune_adt$IDAT <- rownames(immune_adt)

ped <- left_join(cef_ped_df[,c("Epi","Fib","IC","IDAT")], immune_ped, by="IDAT")
adt <- left_join(cef_adt_df[,c("Epi","Fib","IC","IDAT")], immune_adt, by="IDAT")

write.csv(ped, file.path("data", "ped_betas_QCDPB_prc_epidish.csv"), row.names = FALSE)
write.csv(adt, file.path("data", "adt_betas_QCDPB_pr_epidish.csv"), row.names = FALSE)

# ========================
# HEATMAP - PRIMARY REFERENCE 
# ========================
ped <- read.csv(file.path("data", "ped_betas_QCDPB_prc_epidish.csv"))

ss <- read_excel(file.path("ss", PED_META)) %>% 
  filter(Lymph_Node == "F", Batch == "REF")
dim(ss)

ped <- ped %>% filter(IDAT %in% ss$IDAT)

three_df <- ped %>%
  column_to_rownames("IDAT") %>%
  select(Epi, Fib, IC) %>%
  as.matrix()

leuko_df <- ped %>%
  column_to_rownames("IDAT") %>%
  select(NK, CD4T, CD8T, Mono, Neutro, Eosino, B) %>%
  as.matrix()

annot <- ss %>%
  select(IDAT, Clinical_Invasiveness, CC_Cluster, Driver_Group) %>%
  column_to_rownames("IDAT")

annot_cols <- list(
  Clinical_Invasiveness = invasiveness_colors,
  CC_Cluster = cluster_colors,
  Driver_Group = driver_colors
)

p1 <- pheatmap(t(three_df), 
               annotation_col = annot,
               annotation_colors = annot_cols,
               scale = "none",
               clustering_distance_cols = "euclidean",
               show_colnames = FALSE,
               silent = TRUE) 

p2 <- pheatmap(t(leuko_df), 
               annotation_col = annot,
               annotation_colors = annot_cols,
               scale = "none",
               clustering_distance_cols = "euclidean",
               show_colnames = FALSE,
               silent = TRUE)

pdf(file.path("figures", "epidish_combined_heatmap.pdf"), width = 10, height = 6, onefile = FALSE)
grid.arrange(p1[[4]], p2[[4]], 
             nrow = 2, 
             heights = c(0.4, 0.6))
dev.off()

