# FIG 3 - JOINT ANALYSIS 
# Clustering 
# ========================
# LOAD DATA
# ========================
ss <- read_excel(file.path("ss", ADT_META))
primary <- ss %>% dplyr::filter(Lymph_Node == "F")
pca <- readRDS(file.path("data", "adt_betas_QCDPB_pr_30000_pca.rds"))

# ========================
# OPTIMIZE PC NUMBER
# ========================
var <- pca$sdev^2 / sum(pca$sdev^2)
cum_var <- cumsum(var)
num_pcs <- which(cum_var >= 0.90)[1]
cat("90% variance threshold:", num_pcs, "PCs\n")

pcs <- pca$x[, 1:num_pcs]
write.csv(pcs, file.path("data", "adt_unmasked_pca_coords.csv"))

# ========================
# FINAL t-SNE
# ========================
pcs <- read.csv(file.path("data", "adt_unmasked_pca_coords.csv"), row.names = 'X')

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

ss_filtered <- ss %>% filter(IDAT %in% rownames(pcs))
primary_filtered <- primary %>% filter(IDAT %in% rownames(pcs))

cat(sprintf("Running t-SNE on all samples: %d samples\n", nrow(ss_filtered)))
adt_tsne_all <- run_tsne(pcs[ss_filtered$IDAT, ], "adt_tsne_all")

cat(sprintf("Running t-SNE on primary samples: %d samples\n", nrow(primary_filtered)))
adt_tsne_primary <- run_tsne(pcs[primary_filtered$IDAT, ], "adt_tsne_primary")
