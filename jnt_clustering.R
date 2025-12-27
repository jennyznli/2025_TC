# FIG 3 - JOINT ANALYSIS 
# Clustering 
# ========================
# LOAD DATA
# ========================
ss <- read_excel(file.path("ss", JNT_META))
ss_ref <- ss %>% dplyr::filter(Reference == 1) #600
ss_primary_ref <- ss %>% dplyr::filter(Lymph_Node == "F", Reference == 1) #536
ss_primary_all <- ss %>% dplyr::filter(Lymph_Node == "F") #621

pca <- readRDS(file.path("data", "jnt_betas_QCDPB_pr_30000_pca.rds"))

# ========================
# OPTIMIZE PC NUMBER
# ========================
var <- pca$sdev^2 / sum(pca$sdev^2)
cum_var <- cumsum(var)
num_pcs <- which(cum_var >= 0.90)[1]
cat("90% variance threshold:", num_pcs, "PCs\n")

pcs <- as.data.frame(pca$x[, 1:num_pcs])

write.csv(pcs, file.path("data", "jnt_unmasked_pca_coords.csv"))
# ========================
# FINAL t-SNE
# ========================
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

run_tsne(pcs[ss$IDAT,], "jnt_tsne_all")
run_tsne(pcs[ss_ref$IDAT,], "jnt_tsne_ref")
run_tsne(pcs[ss_primary_ref$IDAT,], "jnt_tsne_primary_ref")
run_tsne(pcs[ss_primary_all$IDAT,], "jnt_tsne_primary_all")




