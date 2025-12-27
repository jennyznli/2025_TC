# FIG 1 
# pca on 30k most variable
# generates tsne coordinates for each subset 
# performs consensus clustering 

# ========================
# LOAD DATA
# ========================
source("config.R")

ss <- read_excel(file.path("ss", PED_META))

# primary samples only  
ss_primary_1 <- ss %>% dplyr::filter(Batch == "REF", Primary_Include_In_Analysis == 1)
ss_primary_12 <- ss %>% dplyr::filter(Primary_Include_In_Analysis == 1)

# primary samples w/ LN metastases 
ss_all_1 <- ss %>% dplyr::filter(Batch == "REF", LN_Include_In_Analysis == 1)
ss_all_12 <- ss %>% dplyr::filter(LN_Include_In_Analysis == 1)

pca <- readRDS(file.path("data", "ped_betas_QCDPB_prc_30000_pca.rds"))

# ========================
# OPTIMIZE PC NUMBER
# ========================
var <- pca$sdev^2 / sum(pca$sdev^2)
cum_var <- cumsum(var)

num_pcs <- which(cum_var >= 0.90)[1]
cat("90% variance threshold:", num_pcs, "PCs\n")
pcs <- pca$x[, 1:num_pcs]

write.csv(pcs, file.path("data", "ped_unmasked_pca_coords.csv"))

# ========================
# FINAL t-SNE
# ========================
pcs <- read.csv(file.path("data", "ped_unmasked_pca_coords.csv"), row.names = 'X')

run_tsne <- function(pca_coords, output_name) {
  n_samples <- nrow(pca_coords)
  perplexity <- floor(sqrt(n_samples))
  cat(sprintf("Running t-SNE on %d samples (perplexity = %d)\n", n_samples, perplexity))
  
  set.seed(123)
  tsne <- Rtsne::Rtsne(
    pca_coords,
    dims = 2,
    check_duplicates = FALSE,
    perplexity = perplexity,
    max_iter = 2000,
    pca = FALSE
  )
  
  df <- as.data.frame(tsne$Y)
  colnames(df) <- c("tSNE1", "tSNE2")
  df$IDAT <- rownames(pca_coords)
  write.csv(df, file.path("data", paste0(output_name, "_coords.csv")), row.names = FALSE)
  
  p <- ggplot(df, aes(x = tSNE1, y = tSNE2)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_minimal() +
    labs(
      title = sprintf("t-SNE plot (perplexity = %d)", perplexity),
      x = "t-SNE 1",
      y = "t-SNE 2"
    )
  print(p)
  
  return(df)
}

tsne_primary_1  <- run_tsne(pcs[ss_primary_1$IDAT, ], "ped_tsne_primary_1")
tsne_primary_12 <- run_tsne(pcs[ss_primary_12$IDAT, ], "ped_tsne_primary_12")

tsne_all_1  <- run_tsne(pcs[ss_all_1$IDAT, ], "ped_tsne_all_1")
tsne_all_12 <- run_tsne(pcs[ss_all_12$IDAT, ], "ped_tsne_all_12")

# ========================
# CONSENSUS CLUSTERING - primary reference
# ========================
ss <- read_excel(file.path("ss", PED_META))
pcs <- t(read.csv(file.path("data", "ped_unmasked_pca_coords.csv"), row.names = 'X'))
pcs <- pcs[, ss_primary_1$IDAT]

set.seed(123)
hc_res <- ConsensusClusterPlus::ConsensusClusterPlus(
  d = pcs,
  maxK = 10,
  reps = 1000,
  pItem = 1,
  pFeature = 0.8,
  clusterAlg = "hc",
  distance = "euclidean",
  seed = 123,
  plot = "pdf",
  writeTable = TRUE,
  title = "ped87_hc_euc_pcs_1000r_1i_0.8f"
)

df <- as.data.frame(hc_res[[3]]$consensusClass)
df$IDAT <- rownames(df)
colnames(df) <- c("CC_Cluster", "IDAT")
df$CC_Cluster <- factor(df$CC_Cluster, 
                        levels = c(1, 2, 3),
                        labels = c("LI", "HI", "LEUKO"))
ss <- dplyr::left_join(ss, df, by = "IDAT")

write_xlsx(ss, file.path("ss", "ped87_unmasked_consensus_clusters.xlsx"))


# ========================
# CONSENSUS CLUSTERING - reference
# ========================
ss <- read_excel(file.path("ss", PED_META))
pcs <- read.csv(file.path("data", "ped_unmasked_pca_coords.csv"), row.names = 'X')
pcs <- t(pcs)[, ss_all_1$IDAT]

set.seed(123)
hc_res <- ConsensusClusterPlus::ConsensusClusterPlus(
  d = pcs,
  maxK = 10,
  reps = 1000,
  pItem = 1.0,
  pFeature = 0.8,
  clusterAlg = "hc",
  distance = "euclidean",
  seed = 123,
  plot = "pdf",
  writeTable = TRUE,
  title = "ped100_hc_euc_pcs_1000r_1i_0.8f"
)

df <- as.data.frame(hc_res[[3]]$consensusClass)
df$IDAT <- rownames(df)
colnames(df) <- c("CC_Cluster", "IDAT")
df$CC_Cluster <- factor(df$CC_Cluster, 
                        levels = c(1, 2, 3),
                        labels = c("LI", "HI", "LEUKO"))
ss <- dplyr::left_join(ss, df, by = "IDAT")

write_xlsx(ss, file.path("ss", "ped100_unmasked_consensus_clusters.xlsx"))

# ========================
# CONSENSUS CLUSTERING - all REF & VAL 
# ========================
ss <- read_excel(file.path("ss", PED_META))
pcs <- read.csv(file.path("data", "ped_unmasked_pca_coords.csv"), row.names = 'X')
pcs <- t(pcs)[, ss_all_12$IDAT]

set.seed(123)
hc_res <- ConsensusClusterPlus::ConsensusClusterPlus(
  d = pcs,
  maxK = 10,
  reps = 1000,
  pItem = 1.0,
  pFeature = 0.8,
  clusterAlg = "hc",
  distance = "euclidean",
  seed = 123,
  plot = "pdf",
  writeTable = TRUE,
  title = "ped184_hc_euc_pcs_1000r_1i_0.8f"
)

df <- as.data.frame(hc_res[[3]]$consensusClass)
df$IDAT <- rownames(df)
colnames(df) <- c("CC_Cluster", "IDAT")
df$CC_Cluster <- factor(df$CC_Cluster, 
                        levels = c(1, 2, 3),
                        labels = c("LI", "HI", "LEUKO"))
ss <- dplyr::left_join(ss, df, by = "IDAT")
write_xlsx(ss, file.path("ss", "ped184_unmasked_consensus_clusters.xlsx"))





