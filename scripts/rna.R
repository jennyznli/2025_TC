# ============================================================
# Performs RNA analysis (differential expression) and GSEA. 
# ============================================================
source("config.R")

# ========================
# LOADING DATA
# ========================
df <- read.csv(file.path("data", "20250113_rna_log2cpmfiltered.csv"), row.names = "GENEID")

ss <- read_excel(file.path("ss", PED_META)) %>% 
  filter(Lymph_Node == "F", Batch == "REF")

# ========================
# DIFF EXPRESSION
# ========================
primary <- df[, ss$Sample_ID]

ss$Invasiveness <- factor(ss$Clinical_Invasiveness, levels = c("Low", "High"))
ss$Sex <- factor(ss$Sex)
ss$IC_EpiDISH <- as.numeric(ss$IC_EpiDISH)

design <- model.matrix(~ Sex + IC_EpiDISH + Invasiveness, data = ss)

fit <- lmFit(primary, design)
fit2 <- eBayes(fit)

deg_all_adjusted <- topTable(
  fit2, 
  coef = "InvasivenessHigh",
  n = Inf, 
  adjust.method = "BH", 
  sort.by = "P"
)
head(deg_all_adjusted)
dim(deg_all_adjusted)
# 22545  

n_sig <- sum(deg_all_adjusted$adj.P.Val < 0.05)
cat("\nNumber of DE genes (FDR < 0.05): ", n_sig, "\n") # 6000

cat("  Upregulated:   ", sum(deg_all_adjusted$adj.P.Val < 0.05 & deg_all_adjusted$logFC > 0), "\n")
cat("  Downregulated: ", sum(deg_all_adjusted$adj.P.Val < 0.05 & deg_all_adjusted$logFC < 0), "\n")
# 3637 
# 2363 

saveRDS(deg_all_adjusted, file.path("data", "deg_primary_covar_all.rds"))

deg_stringent_adjusted <- deg_all_adjusted %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)
cat("  Upregulated:   ", sum(deg_stringent_adjusted$adj.P.Val < 0.05 & deg_stringent_adjusted$logFC > 1), "\n")
cat("  Downregulated: ", sum(deg_stringent_adjusted$adj.P.Val < 0.05 & deg_stringent_adjusted$logFC < -1), "\n")
# 1213 
# 492 

saveRDS(deg_stringent_adjusted, file.path("data", "deg_primary_covar_stringent.rds"))
dim(deg_stringent_adjusted)

deg_sig_adj <- deg_all_adjusted %>%
  filter(adj.P.Val < 0.05)
dim(deg_sig_adj)
# 6000
saveRDS(deg_sig_adj, file.path("data", "deg_primary_covar_sig.rds"))

# ========================
# DIFFERENTIAL METHYLATION - TFBS
# ========================
logfc_threshold <- 1
effect_threshold <- 0.05
fdr_threshold <- 0.05

hypo_enrich <- readRDS(file.path("data", "ped_dm_inv_hypo_tfbs_enrichment_leuko.rds"))
hyper_enrich <- readRDS(file.path("data", "ped_dm_inv_hyper_tfbs_enrichment_leuko.rds"))

tf_hyper <- hyper_enrich %>%
  filter(FDR < fdr_threshold) %>%
  pull(dbname) %>%
  unique()
length(tf_hyper)
# 0 
tf_hypo <- hypo_enrich %>%
  filter(FDR < fdr_threshold) %>%
  pull(dbname) %>%
  unique()
length(tf_hypo)
# 264 

# ========================
# SCATTER PLOT - DIFF METH + DIFF EXPRESSION 
# ========================
res <- readRDS(file.path("data", "ped_dm_invasiveness_leuko_pr_bh.rds"))
sig_probes <- res %>%
  filter(Clinical_InvasivenessHigh_BH < fdr_threshold) %>%
  mutate(Methylation = ifelse(Est_Clinical_InvasivenessHigh > 0, "Hyper", "Hypo"))
dim(sig_probes)
# 133851

# get gene annotations from sesameData for probes
regs <- sesameData_getTxnGRanges("hg38", merge2gene = TRUE)
dm_regs <- as.data.frame(sesameData_annoProbes(sig_probes$Probe_ID, regs, genome="hg38", column="gene_name"))
dm_regs$Probe_ID <- rownames(dm_regs)

# methylation - collapse down to gene level 
gene_meth <- res %>%
  left_join(dm_regs %>% select(Probe_ID, Gene = gene_name), by = "Probe_ID") %>%
  filter(!is.na(Gene), Gene != "") %>%
  group_by(Gene) %>%
  summarize(
    n_probes = n(),
    mean_meth_delta = median(Est_Clinical_InvasivenessHigh),
    min_meth_pval = min(Clinical_InvasivenessHigh_BH),
    .groups = "drop"
  )

# merge with diff expression and label based on opposite patterns
integrated <- gene_meth %>%
  left_join(deg_all_adjusted %>% mutate(Gene = rownames(.)), by = "Gene") %>%
  mutate(
    regulation_pattern = case_when(
      mean_meth_delta < -effect_threshold & logFC > logfc_threshold ~ "Activated",
      mean_meth_delta > effect_threshold & logFC < -logfc_threshold ~ "Silenced",
      TRUE ~ NA_character_
    )
  ) 
dim(integrated)
# 21242

tfbs <- integrated %>% filter(Gene %in% c(tf_hyper, tf_hypo))
dim(tfbs)
# 192 

# label top genes significant in both methylation and expression
genes_label <- tfbs %>%
  filter(!is.na(regulation_pattern), min_meth_pval < fdr_threshold, adj.P.Val < fdr_threshold) %>%
  arrange(desc(abs(logFC))) %>%
  head(10) %>%
  pull(Gene)
length(genes_label)

p <- ggplot(tfbs, aes(x = logFC, y = mean_meth_delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = c(-effect_threshold, effect_threshold), 
             linetype = "dashed", color = "gray70", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
             linetype = "dashed", color = "gray70", alpha = 0.5) +
  geom_point(data = tfbs %>% filter(is.na(regulation_pattern)),
             aes(size = n_probes), color = "gray80", alpha = 0.4) +
  geom_point(data = tfbs %>% filter(!is.na(regulation_pattern)),
             aes(color = regulation_pattern, size = n_probes), alpha = 0.8) +
  scale_color_manual(values = c("Activated" = "#d61525ff", "Silenced" = "#0c66bcff")) +
  scale_size_continuous(range = c(1.5, 6), name = "# Probes") +
  ggrepel::geom_text_repel(
    data = subset(tfbs, Gene %in% genes_label),
    aes(label = Gene),
    size = 3.5, fontface = "bold", box.padding = 0.5,
    max.overlaps = 20, min.segment.length = 0
  ) +
  labs(
    x = "Differential Expression (log2 FC)",
    y = "Average DNA Methylation Beta",
    color = "Epigenetic Regulation"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

pdf(file.path("figures", "tfbs_epigenetic_silencing.pdf"), width = 6, height = 5)
print(p)
dev.off()

# ========================
# HEATMAP - TF DIFF EXPRESSION 
# ========================
invasiveness <- factor(ss$Clinical_Invasiveness, levels = c("Low", "High"))

deg_hyper_expressed <- deg_stringent %>%
  filter(rownames(.) %in% tf_hyper)
deg_hypo_expressed <- deg_stringent %>%
  filter(rownames(.) %in% tf_hypo)

deg_tf_expr <- primary[unique(c(rownames(deg_hypo_expressed), rownames(deg_hyper_expressed))), ]

annotation_col <- data.frame(
  Cluster = as.factor(ss$CC_Cluster),
  Invasiveness = as.factor(ss$Clinical_Invasiveness),
  row.names = colnames(deg_tf_expr)
)

annotation_colors <- list(
  Cluster = cluster_colors,
  Invasiveness = invasiveness_colors
)

pdf(file.path("figures", "diff_exp_meth_tf_heatmap.pdf"), width=7, height=3, onefile=FALSE)
pheatmap(deg_tf_expr,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         color = colorRampPalette(rev(brewer.pal(11, "PiYG")))(50),
         border_color = NA
)
dev.off()

# ========================
# GSEA
# ========================
deg_all <- readRDS(file.path(DATA_DIR, "diff_exp_genes.rds"))

ranked_genes <- deg_all %>%
  mutate(
    ranking_metric = -log10(P.Value) * sign(logFC),
    gene_id = rownames(.)
  ) %>%
  arrange(desc(ranking_metric))

gene_ranks <- setNames(ranked_genes$ranking_metric, ranked_genes$gene_id)
h_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  mutate(collection = "H")

msig_list <- split(h_sets$gene_symbol, h_sets$gs_name)

gsea_results <- fgsea(
  pathways = msig_list,
  stats = gene_ranks,
  scoreType = 'std',
  minSize = 15,
  maxSize = 500,
  nPermSimple = 10000,
  eps = 0
) %>%
  as_tibble() %>%  #
  arrange(padj)
print(dim(gsea_results))
# 50  8

gsea_results_filtered <- gsea_results %>%
  filter(padj < 0.1) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 20) %>%
  mutate(
    pathway = pathway %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ")
  )

pdf(file.path(FIG_DIR, "gsea_enrichment_dot_hallmark.pdf"), width=6, height=4, onefile=FALSE)
p <- ggplot(gsea_results_filtered,
            aes(x = NES,
                y = reorder(pathway, NES),
                size = -log10(padj),
                color = NES)) +
  geom_point() +
  scale_color_gradientn(colors = brewer.pal(9, "PuRd")[3:9]) +
  scale_size_continuous(name = "-log10(padj)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Pathway")
plot(p)
dev.off()
