# ============================================================
# This script generates general stat. figures summarizing 
# probe success rates, tissue types, and global mean methylation 
# for the pediatric cohort. 
# ============================================================
source("../config.R")

# ========================
# LOAD DATA
# ========================
ss <- read_excel(file.path("ss", PED_META))

ss$Clinical_Invasiveness <- as.factor(ss$Clinical_Invasiveness)
ss$Clinical_Confidence <- as.factor(ss$Clinical_Confidence)
ss$Fmodel_Invasiveness_Prediction <- as.factor(ss$Fmodel_Invasiveness_Prediction)
ss$Fmodel_Invasiveness_Confidence <- as.numeric(ss$Fmodel_Invasiveness_Confidence)
ss$Fmodel_Invasiveness_Accuracy <- as.numeric(ss$Fmodel_Invasiveness_Accuracy)
ss$Fmodel_Driver_Prediction <- as.factor(ss$Fmodel_Driver_Prediction)
ss$Fmodel_Driver_Confidence <- as.numeric(ss$Fmodel_Driver_Confidence)
ss$Fmodel_Driver_Accuracy <- as.numeric(ss$Fmodel_Driver_Accuracy)
ss$Lymph_Node <- as.factor(ss$Lymph_Node)
ss$Chronological_Age <- as.numeric(ss$Chronological_Age)
ss$M <- as.factor(ss$M)
ss$N <- as.factor(ss$N)
ss$`T` <- as.factor(ss$`T`)
ss$Tissue_Type <- as.factor(ss$Tissue_Type)
ss$Batch <- as.factor(ss$Batch)
ss$Probe_Success_Rate <- as.numeric(ss$Probe_Success_Rate)
ss$IC_EpiDISH <- as.numeric(ss$IC_EpiDISH)
ss$CC_Cluster <- as.factor(ss$CC_Cluster)

# ========================
# TISSUE TYPE BAR PLOT - BATCH
# ========================
tissue_summary <- ss %>% group_by(Batch, Tissue_Type) %>% summarise(count = n(), .groups = 'drop') %>% group_by(Batch) %>% mutate( total = sum(count), proportion = count / total )
print(tissue_summary)

p_stacked <- ggplot(ss, aes(x = Batch, fill = Tissue_Type)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = tissue_colors) +
  labs(title = "Tissue Type Distribution by Batch",
       x = "Batch", y = "Count", fill = "Tissue Type") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures/tissue_type_batch_stacked.pdf",
       p_stacked, width = 4, height = 4.5, device = "pdf")

# ========================
# PROBE SUCCESS RATE BOX PLOTS - BATCH
# ========================
psr_summary <- ss %>%
  group_by(Batch) %>%
  summarise(
    n = n(),
    median = median(Probe_Success_Rate, na.rm = TRUE),
    Q1 = quantile(Probe_Success_Rate, 0.25, na.rm = TRUE),
    Q3 = quantile(Probe_Success_Rate, 0.75, na.rm = TRUE),
    mean = mean(Probe_Success_Rate, na.rm = TRUE),
    sd = sd(Probe_Success_Rate, na.rm = TRUE),
    .groups = "drop"
  )

print(psr_summary)
print(wilcox.test(Probe_Success_Rate ~ Batch, data = ss))
# Wilcoxon rank sum test with continuity correction
# 
# data:  Probe_Success_Rate by Batch
# W = 7271, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
# 
p_batch <- ggboxplot(ss, x = "Batch", y = "Probe_Success_Rate",
                     fill = "Batch", legend = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
  labs(title = "Probe Success Rate by Batch",
       x = "Batch", y = "Probe Success Rate") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = batch_colors)

ggsave("figures/psr_batch_boxplots.pdf", p_batch, width = 3.5, height = 4.5, device = "pdf")

# ========================
# GLOBAL MEAN BETA BOX PLOTS - CLUSTER, INVASIVENESS, AND DRIVER GROUP
# =======================
library(patchwork)

betas <- readRDS(file.path("data", PED_BETAS))

ss_ref <- ss %>%
  filter(Batch == "REF", Lymph_Node == "F")
betas_ref <- betas[, ss_ref$IDAT]
ss_ref$Global_Means <- colMeans(betas_ref)

p_inv <- ggplot(ss_ref, aes(x = Clinical_Invasiveness, y = Global_Means, fill = Clinical_Invasiveness)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = invasiveness_colors) +
  labs(title = "Invasiveness", x = NULL, y = "Global Mean DNAm") +
  theme_minimal() +
  theme(legend.position = "none")

p_cluster <- ggplot(ss_ref, aes(x = CC_Cluster, y = Global_Means, fill = CC_Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cluster", x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_driver <- ss_ref %>%
  filter(Driver_Group != "Indeterminate") %>%
  ggplot(aes(x = Driver_Group, y = Global_Means, fill = Driver_Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = driver_colors) +
  labs(title = "Driver", x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

combined <- wrap_plots(p_inv, p_cluster, p_driver, nrow = 1, widths = c(1, 1.5, 2))

ggsave("figures/primary_globalmeans_combined.pdf", combined, width = 7, height = 3.5)


# ========================
# PROBE SUCCESS RATE VS. CONFIDENCE OF ALL CLASSIFIERS 
# ========================
make_plot <- function(df, acc_col, conf_col, title = NULL) {
  stopifnot(acc_col %in% colnames(df))
  stopifnot(conf_col %in% colnames(df))
  
  acc_raw <- df[[acc_col]]
  conf_raw <- df[[conf_col]]
  
  if (!is.numeric(conf_raw)) {
    suppressWarnings(conf_num <- as.numeric(conf_raw))
    if (all(is.na(conf_num))) {
      stop(conf_col, " could not be coerced to numeric.")
    }
  } else {
    conf_num <- conf_raw
  }
  
  if (is.logical(acc_raw)) {
    acc_chr <- ifelse(acc_raw, "Incorrect", "Correct")
  } else if (is.numeric(acc_raw)) {
    if (!all(na.omit(acc_raw) %in% c(0,1))) stop(acc_col, " must contain only 0/1.")
    acc_chr <- ifelse(acc_raw == 1, "Incorrect", "Correct")
  } else if (is.character(acc_raw)) {
    acc_chr <- case_when(
      tolower(acc_raw) %in% c("0","correct") ~ "Correct",
      tolower(acc_raw) %in% c("1","incorrect") ~ "Incorrect",
      TRUE ~ NA_character_
    )
    if (all(is.na(acc_chr))) stop(acc_col, " has invalid values. Expected 0/1 or Correct/Incorrect.")
  } else {
    stop(acc_col, " must be logical, numeric, or character.")
  }
  
  df2 <- tibble(acc_bin = acc_chr, conf_val = conf_num) %>%
    filter(!is.na(acc_bin), !is.na(conf_val))
  
  
  wilcox_res <- wilcox.test(conf_val ~ acc_bin, data = df2, alternative = "greater", exact = FALSE)
  
  p_label <- if (wilcox_res$p.value < 0.001) {
    "p < 0.001"
  } else {
    sprintf("p = %.3f", wilcox_res$p.value)
  }
  
  ggplot(df2, aes(x = acc_bin, y = conf_val, fill = acc_bin)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    annotate("text", x = 1.5, y = 0.95, label = p_label, size = 3.5) +
    scale_fill_manual(values = c("Correct"="#4DAF4A","Incorrect"="#E41A1C")) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    labs(x = NULL, y = "Prediction Confidence", title = title) +
    theme_minimal(base_size = 12) +
    theme( legend.position = "none", plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
}

ss_ref <- read_excel(file.path("ss", PED_META)) %>%
  filter(Batch == "REF", Lymph_Node == "F")

p1 <- make_plot(
  df = ss_ref,
  acc_col = "LOOCV_Invasiveness_Accuracy",
  conf_col = "LOOCV_Invasiveness_Confidence",
  title = "Reference – Invasiveness"
)

ss_val <- read_excel(file.path("ss", PED_META)) %>%
  filter(Batch == "VAL", Lymph_Node == "F")

p2 <- make_plot(
  df = ss_val,
  acc_col = "Fmodel_Invasiveness_Accuracy",
  conf_col = "Fmodel_Invasiveness_Confidence",
  title = "Validation – Invasiveness"
)

ss_ref <- ss_ref %>% filter(Driver_Group != "Indeterminate")
ss_val <- ss_val %>% filter(Driver_Group != "Indeterminate")

p3 <- make_plot(
  df = ss_ref,
  acc_col = "LOOCV_Driver_Accuracy",
  conf_col = "LOOCV_Driver_Confidence",
  title = "Reference – Driver"
)

p4 <- make_plot(
  df = ss_val,
  acc_col = "Fmodel_Driver_Accuracy",
  conf_col = "Fmodel_Driver_Confidence",
  title = "Validation – Driver"
)

final_plot <- p1 | p2 | p3 | p4

ggsave(filename = file.path("figures", "all_accuracy_vs_confidence.pdf"), plot = final_plot, width = 10, height = 4)



