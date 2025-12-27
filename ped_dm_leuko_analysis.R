source(file.path("config.R"))

# ========================
# BH CORRECTION
# ========================
### INVASIVENESS 
res <- readRDS(here("data", "ped_dm_invasiveness_leuko.rds"))
res$Clinical_InvasivenessHigh_BH <- p.adjust(res$Pval_Clinical_InvasivenessHigh, method = "BH")
saveRDS(res, file.path("data", "ped_dm_invasiveness_leuko_bh.rds"))

### CLUSTERS
res1 <- readRDS(here("data", "ped_dm_cluster_li.rds"))
res2 <- readRDS(here("data", "ped_dm_cluster_leuko.rds"))

res1$Pval_CC_ClusterHI_BH <- p.adjust(res1$Pval_CC_ClusterHI, method = "BH")
res1$Pval_CC_ClusterLEUKO_BH <- p.adjust(res1$Pval_CC_ClusterLEUKO, method = "BH")
res2$Pval_CC_ClusterHI_BH <- p.adjust(res2$Pval_CC_ClusterHI, method = "BH")
res2$Pval_CC_ClusterLI_BH <- p.adjust(res2$Pval_CC_ClusterLI, method = "BH")

saveRDS(res1, file.path("data", "ped_dm_cluster_li_bh.rds"))
saveRDS(res2, file.path("data", "ped_dm_cluster_leuko_bh.rds"))

# ========================
# 2A. INVASIVENESS VOLCANO
# ========================
res <- readRDS(here("data", "ped_dm_invasiveness_leuko_bh.rds"))
dim(res)
dim(res %>% filter(Clinical_InvasivenessHigh_BH < 0.05))

res <- res %>%
  mutate(
    threshold = case_when(
      Clinical_InvasivenessHigh_BH < 0.05 & Est_Clinical_InvasivenessHigh >  0.2  ~ "Hyper",
      Clinical_InvasivenessHigh_BH < 0.05 & Est_Clinical_InvasivenessHigh < -0.2  ~ "Hypo",
      TRUE ~ "NS"
    )
  )
table(res$threshold)
# Hyper   Hypo     NS 
# 150   4706 897196 

sig_points <- subset(res, threshold != "NS")
nonsig_points <- subset(res, threshold == "NS")

saveRDS(res, file.path("data", "ped_dm_invasiveness_leuko_bh_anno.rds"))

p <- ggplot() +
  geom_bin2d(data = nonsig_points,
             aes(x = Est_Clinical_InvasivenessHigh, y = -log10(Clinical_InvasivenessHigh_BH)),
             bins = 250) +
  scale_fill_gradient(low = "lightgray",
                      high = "black",
                      trans = "log10") + # log transform to better show density differences
  geom_point(data = sig_points,
             aes(x = Est_Clinical_InvasivenessHigh, y = -log10(Clinical_InvasivenessHigh_BH),
                 color = threshold),
             size = 0.1,
             alpha = 0.3) +
  scale_color_manual(values = c(
    "Hypo" = "#0c66bcff",
    "Hyper" = "#d61525ff"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "#575757ff") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "#575757ff") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#575757ff") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey90"),
    legend.position = "right",
    legend.title = element_blank()
  ) +
  labs(
    x = "Methylation Difference",
    y = "-log10(P-value)"
  ) +
  guides(fill = "none")

pdf(file.path("figures", "ped_invasiveness_volcano_leuko.pdf"), width=5, height=5, onefile=FALSE)
plot(p)
dev.off()

# ========================
# 2B. HEATMAP 88
# ========================
betas <- readRDS(file.path("data", PED_BETAS))
ss_primary <- read_excel(file.path("ss", PED_META)) %>% 
  filter(Lymph_Node == "F", Primary_Include_In_Analysis == "1", Batch == "REF")
res <- readRDS(file.path("data", "ped_dm_invasiveness_leuko_bh.rds"))

# Filter for significant probes
sig_probes <- res %>%
  dplyr::filter(Clinical_InvasivenessHigh_BH < 0.05,
                abs(Est_Clinical_InvasivenessHigh) > 0.2) %>%
  mutate(Methylation = ifelse(Est_Clinical_InvasivenessHigh > 0.2, "Hyper", "Hypo"))

# Create annotations
annotation_row <- data.frame(
  Methylation = factor(sig_probes$Methylation),
  row.names = sig_probes$Probe_ID
)

annotation_col <- data.frame(
  Clinical_Invasiveness = factor(ss_primary$Clinical_Invasiveness),
  CC_Cluster = factor(ss_primary$CC_Cluster),
  row.names = ss_primary$Sample_ID
) %>% arrange(Clinical_Invasiveness)

betas2 <- betas[sig_probes$Probe_ID, ss_primary$IDAT]
dim(betas2)
colnames(betas2) <- ss_primary$Sample_ID[match(colnames(betas2), ss_primary$IDAT)]

annotation_colors <- list(
  Methylation = methylation_colors,
  Clinical_Invasiveness = invasiveness_colors,
  CC_Cluster = cluster_colors
)

pdf(file.path("figures", "ped_dm_invasiveness_heatmap_leuko.pdf"), width = 13, height = 12, onefile = FALSE)
pheatmap(betas2,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#0070ffff", "white", "#f70000ff"))(50),
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_legend = FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         legend = FALSE)
dev.off()

# ========================
# S2A. HEATMAP 101
# ========================
betas <- readRDS(file.path("data", PED_BETAS))

ss <- read_excel(file.path("ss", PED_META)) %>% 
  filter(LN_Include_In_Analysis == "1", Batch == "REF")
res <- readRDS(file.path("data", "ped_dm_invasiveness_leuko_bh.rds"))

sig_probes <- res %>%
  dplyr::filter(Clinical_InvasivenessHigh_BH < 0.05,
                abs(Est_Clinical_InvasivenessHigh) > 0.2) %>%
  mutate(Methylation = ifelse(Est_Clinical_InvasivenessHigh > 0.2, "Hyper", "Hypo"))

annotation_col <- data.frame(
  Clinical_Invasiveness = factor(ss$Clinical_Invasiveness),
  Cluster_Group = factor(ss$CC_Cluster),
  Lymph_Node = factor(ss$Lymph_Node),
  row.names = ss$Sample_ID 
) %>% arrange(Clinical_Invasiveness)

betas2 <- betas[sig_probes$Probe_ID, ss$IDAT]
dim(betas2)
colnames(betas2) <- ss$Sample_ID[match(colnames(betas2), ss$IDAT)]

# Create annotations
annotation_row <- data.frame(
  Methylation = factor(sig_probes$Methylation),
  row.names = sig_probes$Probe_ID
)

annotation_colors <- list(
  Cluster_Group = cluster_colors,
  Clinical_Invasiveness = invasiveness_colors,
  Lymph_Node = lymph_node_colors,
  Methylation = methylation_colors
)

pdf(file.path("figures", "ped_ln_dm_invasiveness_heatmap_leuko.pdf"), width = 15, height = 12, onefile = FALSE)
pheatmap(betas2,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#0070ffff", "white", "#f70000ff"))(50),
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_legend = FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         legend = FALSE)
dev.off()

# ========================
# 2D-E, S2E-G. ENRICHMENT PLOTS - INVASIVENESS
# ========================
hyper <- sig_probes %>%
  filter(Methylation == "Hyper")
hypo <- sig_probes %>%
  filter(Methylation == "Hypo")

# LOW vs HIGH — TFBS HYPO
x <- testEnrichment(hypo$Probe_ID, "TFBSconsensus", platform = "EPIC", universe = res$Probe_ID)
pdf(file.path("figures", "ped_dm_inv_hypo_tfbs_enrichment_leuko.pdf"),
    width = 3, height = 4, onefile = FALSE)
plotDotFDR(x, n_min = 30, n_max = 30)
dev.off()
saveRDS(x, file.path("data", "ped_dm_inv_hypo_tfbs_enrichment_leuko.rds"))

# LOW vs HIGH — TFBS HYPER
y <- testEnrichment(hyper$Probe_ID, "TFBSconsensus", platform = "EPIC", universe = res$Probe_ID)
pdf(file.path("figures", "ped_dm_inv_hyper_tfbs_enrichment_leuko.pdf"),
    width = 3, height = 4, onefile = FALSE)
plotDotFDR(y, n_min = 30, n_max = 30)
dev.off()
saveRDS(y, file.path("data", "ped_dm_inv_hyper_tfbs_enrichment_leuko.rds"))

# ========================
# S2B. TOTAL DMP BARPLOT
# ========================
res1 <- readRDS(file.path("data", "ped_dm_cluster_li_bh.rds"))
res2 <- readRDS(file.path("data", "ped_dm_cluster_leuko_bh.rds"))

# LI vs HI (HI in res1)
cat(dim(filter(res1, Pval_CC_ClusterHI_BH < 0.05, Est_CC_ClusterHI > 0.2))[1], "hyper (LI vs HI)\n")
cat(dim(filter(res1, Pval_CC_ClusterHI_BH < 0.05, Est_CC_ClusterHI < -0.2))[1], "hypo (LI vs HI)\n\n")
# 2502 hyper (LI vs HI)
# 16623 hypo (LI vs HI)

# LI vs LEUKO (HIL in res1)
cat(dim(filter(res1, Pval_CC_ClusterLEUKO_BH < 0.05, Est_CC_ClusterLEUKO > 0.2))[1], "hyper (LI vs LEUKO)\n")
cat(dim(filter(res1, Pval_CC_ClusterLEUKO_BH < 0.05, Est_CC_ClusterLEUKO < -0.2))[1], "hypo (LI vs LEUKO)\n\n")
# 46201 hyper (LI vs LEUKO)
# 37525 hypo (LI vs LEUKO)

# LEUKO vs HI (HI in res2)
cat(dim(filter(res2, Pval_CC_ClusterHI_BH < 0.05, Est_CC_ClusterHI > 0.2))[1], "hyper (LEUKO vs HI)\n")
cat(dim(filter(res2, Pval_CC_ClusterHI_BH < 0.05, Est_CC_ClusterHI < -0.2))[1], "hypo (LEUKO vs HI)\n")
# 16177 hyper (LEUKO vs HI)
# 44113 hypo (LEUKO vs HI)

sig_results <- list(
  "LI / LEUKO" = get_sig_probes(res1, "Pval_CC_ClusterLEUKO_BH", "Est_CC_ClusterLEUKO"),
  "LI / HI" = get_sig_probes(res1, "Pval_CC_ClusterHI_BH", "Est_CC_ClusterHI"),
  "LEUKO / HI" = get_sig_probes(res2, "Pval_CC_ClusterHI_BH", "Est_CC_ClusterHI")
)

data <- bind_rows(
  lapply(names(sig_results), function(group) {
    data.frame(
      Group = group,
      table(sig_results[[group]]$Methylation)
    ) %>%
      arrange(desc(Var1))
  })
)
total_dmps <- data.frame(
  Group = names(sig_results),
  Total = sapply(sig_results, nrow),
  stringsAsFactors = FALSE
) %>%
  mutate(Label = format(Total, big.mark = ","))

data$Methylation <- factor(data$Var1, levels = c("Hyper", "Hypo"))

pdf(file.path("figures", "ped_dm_cluster_composition_leuko.pdf"), width=4, height=4, onefile=FALSE)
p <- ggplot(data, aes(x = Group, y = Freq/10000)) +
  geom_bar(aes(fill = Var1), stat = "identity", position = "stack") +
  scale_fill_manual(values = methylation_colors) +
  geom_text(data = total_dmps,
            aes(y = Total/10000, label = Label),
            vjust = -0.3,
            size = 3.5) +
  scale_y_continuous(name = expression("DMPs [" * 10^4 * "]"),
                     limits = c(0, max(total_dmps$Total/10000) * 1.1),
                     breaks = seq(0, ceiling(max(total_dmps$Total/10000)), by = 2)) +
  theme_minimal() +
  theme(legend.title = element_text(size = 10),
        panel.grid.major.x = element_blank()) +
  labs(fill = "Methylation",
       x = "Cluster Comparison",
       y = expression("DMPs [" * 10^4 * "]"))
plot(p)
dev.off()

# ========================
# S2C. GLOBAL MEANS CLUSTERS
# ========================
betas <- readRDS(file.path("data", PED_BETAS))
ss <- read_excel(file.path("ss", PED_META))
ss_primary <- ss %>% filter(Primary_Include_In_Analysis == 1, Batch == "REF")
betas1 <- betas[, ss_primary$IDAT]

ss_primary$Global_Means = colMeans(betas1)

wcc <- pairwise.wilcox.test(ss_primary$Global_Means, ss_primary$CC_Cluster,
                            p.adjust.method = "BH")
print(wcc)
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  ss_primary$Global_Means and ss_primary$CC_Cluster 
# 
# HI      LEUKO  
# LEUKO 0.00021 -      
#   LI    0.00044 0.58988
# 
# P value adjustment method: BH 

pdf(file.path("figures", "ped_cluster_global_means_wc_leuko.pdf"), width=4.5, height=4.5, onefile=FALSE)
p <- ggplot(ss_primary, aes(x = CC_Cluster, y = Global_Means, fill = CC_Cluster)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("HI", "LEUKO"),
                                        c("HI", "LI")),
                     label = "p") +
  scale_fill_manual(values = cluster_colors) +
  labs(x = "Cluster",
       y = "Global Mean Methylation",
       fill = "Cluster") +
  theme_minimal()
plot(p)
dev.off()
#
# ========================
# S2D. GLOBAL MEANS DRIVERS
# ========================
# Filter out Indeterminate category
ss_filtered <- ss_primary %>%
  filter(Driver_Group != "Indeterminate")

wcc <- pairwise.wilcox.test(ss_filtered$Global_Means, ss_filtered$Driver_Group,
                            p.adjust.method = "BH")
print(wcc)

# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  ss_filtered$Global_Means and ss_filtered$Driver_Group 
# 
# BRAF V600E DICER1  Kinase Fusion
# DICER1        0.87229    -       -            
#   Kinase Fusion 0.87229    0.87229 -            
#   Ras-like      0.00042    0.08569 0.00042      
# 
# P value adjustment method: BH 

pdf(file.path("figures", "ped_driver_global_means_wc_leuko.pdf"), width=5.5, height=5, onefile=FALSE)
p <- ggplot(ss_filtered, aes(x = Driver_Group, y = Global_Means, fill = Driver_Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("BRAF V600E", "Ras-like"),
                                        c("Kinase Fusion", "Ras-like"),
                                        c("DICER1", "Ras-like")),
                     label = "p",
                     hide.ns = TRUE,
                     p.cutoff = 0.01) +
  scale_fill_manual(values = driver_colors) +
  labs(x = "Cluster",
       y = "Global Mean Methylation",
       fill = "Cluster") +
  theme_minimal()
plot(p)
dev.off()

# ========================
# 2D. ENRICHMENT - CLUSTERS, TISSUES
# ========================
res1 <- readRDS(file.path("data", "ped_dm_cluster_li_bh.rds"))
res2 <- readRDS(file.path("data", "ped_dm_cluster_leuko_bh.rds"))

# LEUKO vs HI — TISSUE HYPER
c <- testEnrichment(
  res2$Probe_ID[res2$Pval_CC_ClusterHI_BH < 0.05 & res2$Est_CC_ClusterHI > 0.2],
  "tissueSignature", platform = "EPIC",
  universe = res2$Probe_ID
)
pdf(file.path("figures", "hil_vs_hi_hyper_tissue_enrichment_leuko.pdf"), width = 4.75, height = 2, onefile = FALSE)
plotDotFDR(c, n_max = 10, n_min = 10)
dev.off()

# LEUKO vs HI — TISSUE HYPO
d <- testEnrichment(
  res2$Probe_ID[res2$Pval_CC_ClusterHI_BH < 0.05 & res2$Est_CC_ClusterHI < -0.2],
  "tissueSignature", platform = "EPIC",
  universe = res2$Probe_ID
)
pdf(file.path("figures", "hil_vs_hi_hypo_tissue_enrichment_leuko.pdf"), width = 5, height = 2, onefile = FALSE)
plotDotFDR(d, n_max = 10, n_min = 10)
dev.off()

# LEUKO vs LI — TISSUE HYPER
a <- testEnrichment(
  res2$Probe_ID[res2$Pval_CC_ClusterLI_BH < 0.05 & res2$Est_CC_ClusterLI > 0.2],
  "tissueSignature", platform = "EPIC",
  universe = res2$Probe_ID
)
pdf(file.path("figures", "hil_vs_li_hyper_tissue_enrichment_leuko.pdf"), width = 4.75, height = 2, onefile = FALSE)
plotDotFDR(a, n_max = 10, n_min = 10)
dev.off()

# LEUKO vs LI — TISSUE HYPO
b <- testEnrichment(
  res2$Probe_ID[res2$Pval_CC_ClusterLI_BH < 0.05 & res2$Est_CC_ClusterLI < -0.2],
  "tissueSignature", platform = "EPIC",
  universe = res2$Probe_ID
)
pdf(file.path("figures", "hil_vs_li_hypo_tissue_enrichment_leuko.pdf"), width = 4.75, height = 2, onefile = FALSE)
plotDotFDR(b, n_max = 10, n_min = 10)
dev.off()

# ========================
# S2F,G. ENRICHMENT - CLUSTERS, TFBS
# ========================
# LEUKO vs LI — TFBS HYPER
a <- testEnrichment(
  res2$Probe_ID[res2$Pval_CC_ClusterLI_BH < 0.05 & res2$Est_CC_ClusterLI > 0.2],
  "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path("figures", "hil_vs_li_hyper_tfbs_enrichment_leuko.pdf"), width = 2.9, height = 3.9, onefile = FALSE)
plotDotFDR(a, n_min = 30, n_max = 30)
dev.off()

# LEUKO vs HI — TFBS HYPER
b <- testEnrichment(
  res2$Probe_ID[res2$Pval_CC_ClusterHI_BH < 0.05 & res2$Est_CC_ClusterHI > 0.2], 
  "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path("figures", "hil_vs_hi_hyper_tfbs_enrichment_leuko.pdf"), width = 2.8, height = 4, onefile = FALSE)
plotDotFDR(b, n_min = 30, n_max = 30)
dev.off()

# LEUKO vs LI — TFBS HYPO
c <- testEnrichment(
  res2$Probe_ID[res2$Pval_CC_ClusterLI_BH < 0.05 & res2$Est_CC_ClusterLI < -0.2],
  "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path("figures", "hil_vs_li_hypo_tfbs_enrichment_leuko.pdf"), width = 2.9, height = 4, onefile = FALSE)
plotDotFDR(c, n_min = 30, n_max = 30)
dev.off()

# LEUKO vs HI — TFBS HYPO
d <- testEnrichment(
  res2$Probe_ID[res2$Pval_CC_ClusterHI_BH < 0.05 & res2$Est_CC_ClusterHI < -0.2],
  "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path("figures", "hil_vs_hi_hypo_tfbs_enrichment_leuko.pdf"), width = 2.8, height = 4, onefile = FALSE)
plotDotFDR(d, n_min = 30, n_max = 30)
dev.off()


# LEUKO vs LI — TFBS HYPO
c <- testEnrichment(
  res2$Probe_ID[res2$Pval_CC_ClusterLI_BH < 0.05 & res2$Est_CC_ClusterLI < -0.2],
  "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path("figures", "hil_vs_li_hypo_tfbs_enrichment_leuko.pdf"), width = 2.9, height = 4, onefile = FALSE)
plotDotFDR(c, n_min = 30, n_max = 30)
dev.off()

# ========================
# HI VS LI 
# ========================

# HI vs LI — TISSUE HYPO
a <- testEnrichment(
  res1$Probe_ID[res1$Pval_CC_ClusterHI_BH < 0.05 & res1$Est_CC_ClusterHI < -0.2],
  "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path("figures", "li_vs_hi_hypo_tissue_enrichment_leuko.pdf"), width = 2.8, height = 4, onefile = FALSE)
plotDotFDR(a, n_min = 30, n_max = 30)
dev.off()

# LEUKO vs LI — TISSUE HYPER
b <- testEnrichment(
  res1$Probe_ID[res1$Pval_CC_ClusterLI_BH < 0.05 & res1$Est_CC_ClusterLI > 0.2],
  "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path("figures", "li_vs_hi_hyper_tissue_enrichment_leuko.pdf"), width = 2.9, height = 4, onefile = FALSE)
plotDotFDR(b, n_min = 30, n_max = 30)
dev.off()

# HI vs LI — TFBS HYPO
d <- testEnrichment(
  res1$Probe_ID[res1$Pval_CC_ClusterHI_BH < 0.05 & res1$Est_CC_ClusterHI < -0.2],
  "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path("figures", "li_vs_hi_hypo_tfbs_enrichment_leuko.pdf"), width = 2.8, height = 4, onefile = FALSE)
plotDotFDR(d, n_min = 30, n_max = 30)
dev.off()

# LEUKO vs LI — TFBS HYPER
c <- testEnrichment(
  res1$Probe_ID[res1$Pval_CC_ClusterLI_BH < 0.05 & res1$Est_CC_ClusterLI > 0.2],
  "TFBS", platform = "EPIC", universe = res2$Probe_ID
)
pdf(file.path("figures", "li_vs_hi_hyper_tfbs_enrichment_leuko.pdf"), width = 2.9, height = 4, onefile = FALSE)
plotDotFDR(c, n_min = 30, n_max = 30)
dev.off()

