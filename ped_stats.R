# ============================================================
# GENERAL STATS - SUPPLEMENTARY FIG 1, 4
# ============================================================
source("config.R")

meta_all <- read_excel(file.path("ss", PED_META))

# ========================
# LOAD DATA
# ========================
ss <- read_excel(file.path("ss", PED_META)) %>% 
  filter(LN_Include_In_Analysis == 1)

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
# FIG S1A. PROBE SUCCESS RATE BY BATCH
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

ggsave("figures/psr_batch_boxplots.pdf",
       p_batch, width = 3.5, height = 4.5, device = "pdf")


# ========================
# FIG S1B. TISSUE TYPE DISTRIBUTION BY BATCH
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





betas <- readRDS(file.path("data", PED_BETAS))

ss_ref <- ss %>%
  filter(Batch == "REF", Primary_Include_In_Analysis == 1)

betas_ref <- betas[, ss_ref$IDAT]
ss_ref$Global_Means <- colMeans(betas_ref)

long_df <- bind_rows(
  ss_ref %>%
    mutate(Grouping = "Invasiveness",
           Group = Clinical_Invasiveness),
  
  ss_ref %>%
    mutate(Grouping = "Cluster",
           Group = CC_Cluster),
  
  ss_ref %>%
    filter(Driver_Group != "Indeterminate") %>%
    mutate(Grouping = "Driver",
           Group = Driver_Group)
) %>%
  filter(!is.na(Group))
p_all <- ggplot(long_df,
                aes(x = Group, y = Global_Means, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_wrap(~Grouping, scales = "free_x") +
  scale_fill_manual(values = c(
    invasiveness_colors,
    cluster_colors,
    driver_colors
  )) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(y = "Global Mean DNA Methylation",
       x = NULL)

ggsave(file.path("figures", "all_primary_globalmeans.pdf"),
       p_all, width = 6, height = 3)








library(patchwork)

p_inv <- ggplot(
  ss_ref,
  aes(x = Clinical_Invasiveness, y = Global_Means, fill = Clinical_Invasiveness)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = invasiveness_colors) +
  labs(title = "Invasiveness", x = NULL, y = "Global Mean DNAm") +
  theme_minimal() +
  theme(legend.position = "none")

p_cluster <- ggplot(
  ss_ref,
  aes(x = CC_Cluster, y = Global_Means, fill = CC_Cluster)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cluster", x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_driver <- ss_ref %>%
  filter(Driver_Group != "Indeterminate") %>%
  ggplot(aes(x = Driver_Group, y = Global_Means, fill = Driver_Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = driver_colors) +
  labs(title = "Driver", x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

library(patchwork)

combined <- wrap_plots(
  p_inv, p_cluster, p_driver,
  nrow = 1,
  widths = c(1, 1.5, 2)  # relative widths
)

ggsave("figures/primary_globalmeans_combined.pdf", combined, width = 7, height = 3.5)


# ========================
# GLOBAL MEANS - INVASIVENESS
# ========================
betas <- readRDS(file.path("data", PED_BETAS))

ss_primary <- ss %>%
  filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>%
  mutate(CC_Cluster = droplevels(CC_Cluster))

betas1 <- betas[, ss_primary$IDAT]
ss_primary$Global_Means <- colMeans(betas1)

sigC <- get_sig_comparisons(ss_primary, "Clinical_Invasiveness", "Global_Means")

p_S2C <- plot_group_box(
  df = ss_primary,
  group_var = "Clinical_Invasiveness",
  value_var = "Global_Means",
  fill_vals = invasiveness_colors,
  sig_comps = sigC
)

ggsave(file.path("figures", "primary_globalmeans_inv_wc.pdf"),
       p_S2C, width = 3, height = 4.5)

ss_primary %>%
  dplyr::group_by(Clinical_Invasiveness) %>%
  dplyr::summarize(
    n = sum(!is.na(Global_Means)),
    median = median(Global_Means, na.rm = TRUE),
    q1 = quantile(Global_Means, 0.25, na.rm = TRUE),
    q3 = quantile(Global_Means, 0.75, na.rm = TRUE),
    IQR = IQR(Global_Means, na.rm = TRUE),
    min = min(Global_Means, na.rm = TRUE),
    max = max(Global_Means, na.rm = TRUE),
    .groups = "drop"
  )

# ========================
# GLOBAL MEANS - CLUSTER
# ========================
betas <- readRDS(file.path("data", PED_BETAS))

ss_primary <- ss %>%
  filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>%
  mutate(CC_Cluster = droplevels(CC_Cluster))

betas1 <- betas[, ss_primary$IDAT]
ss_primary$Global_Means <- colMeans(betas1)

sigC <- get_sig_comparisons(ss_primary, "CC_Cluster", "Global_Means")

p_S2C <- plot_group_box(
  df = ss_primary,
  group_var = "CC_Cluster",
  value_var = "Global_Means",
  fill_vals = cluster_colors,
  sig_comps = sigC
)

ggsave(file.path("figures", "primary_globalmeans_cluster_wc.pdf"),
       p_S2C, width = 4, height = 4.5)

ss_primary %>%
  dplyr::group_by(CC_Cluster) %>%
  dplyr::summarize(
    n = sum(!is.na(Global_Means)),
    median = median(Global_Means, na.rm = TRUE),
    q1 = quantile(Global_Means, 0.25, na.rm = TRUE),
    q3 = quantile(Global_Means, 0.75, na.rm = TRUE),
    IQR = IQR(Global_Means, na.rm = TRUE),
    min = min(Global_Means, na.rm = TRUE),
    max = max(Global_Means, na.rm = TRUE),
    .groups = "drop"
  )


# ========================
# FIG S2D. GLOBAL MEANS BY DRIVER GROUP
# ========================
ss_driver <- ss %>%
  filter(
    Batch == "REF",
    Primary_Include_In_Analysis == 1,
    Driver_Group != "Indeterminate"
  ) %>%
  mutate(Driver_Group = factor(Driver_Group))

betas2 <- betas[, ss_driver$IDAT]
ss_driver$Global_Means <- colMeans(betas2)

sigD <- get_sig_comparisons(ss_driver, "Driver_Group", "Global_Means")

p_S2D <- plot_group_box(
  df = ss_driver,
  group_var = "Driver_Group",
  value_var = "Global_Means",
  fill_vals = driver_colors,
  sig_comps = sigD
)

ggsave(file.path("figures", "primary_globalmeans_driver_wc.pdf"),
       p_S2D, width = 4, height = 4.5)

ss_primary %>%
  dplyr::group_by(Driver_Group) %>%
  dplyr::summarize(
    n = sum(!is.na(Global_Means)),
    median = median(Global_Means, na.rm = TRUE),
    q1 = quantile(Global_Means, 0.25, na.rm = TRUE),
    q3 = quantile(Global_Means, 0.75, na.rm = TRUE),
    IQR = IQR(Global_Means, na.rm = TRUE),
    min = min(Global_Means, na.rm = TRUE),
    max = max(Global_Means, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# STATS FOR PAPER
# ============================================================
ss <- read_excel(file.path("ss", PED_META)) %>% 
  filter(LN_Include_In_Analysis == 1)

ss$Clinical_Invasiveness <- as.factor(ss$Clinical_Invasiveness)
ss$Clinical_Confidence <- as.factor(ss$Clinical_Confidence)
ss$Fmodel_Invasiveness_Predicted <- as.factor(ss$Fmodel_Invasiveness_Predicted)
ss$Fmodel_Invasiveness_Confidence <- as.numeric(ss$Fmodel_Invasiveness_Confidence)
ss$Fmodel_Invasiveness_Accuracy <- as.numeric(ss$Fmodel_Invasiveness_Accuracy)
ss$Fmodel_Driver_Predicted <- as.factor(ss$Fmodel_Driver_Predicted)
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

length(unique(ss$Sample_ID)) #184
ss$Patient_ID <- substr(ss$Sample_ID, 1, 7)
length(unique(ss$Patient_ID)) #169

ss %>% filter(LN_Include_In_Analysis == 1) %>% select(Batch) %>% table()

ss %>% filter(LN_Include_In_Analysis == 1, Batch == "REF") %>% select(Tissue_Type) %>% table()
ss %>% filter(LN_Include_In_Analysis == 1, Batch == "VAL") %>% select(Tissue_Type) %>% table()

ss %>% filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>% select(Sex) %>% table()
ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>% select(Sex) %>% table()
ss %>% filter(Batch == "REF", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(Sex) %>% table()
ss %>% filter(Batch == "VAL", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(Sex) %>% table()


ss %>% filter(Primary_Include_In_Analysis == 1, Batch == "REF") %>%
  summarize(
    mean = mean(Chronological_Age),
    sd = sd(Chronological_Age, na.rm = TRUE),
    min = min(Chronological_Age),
    max = max(Chronological_Age)
  )
ss %>%
  filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>%
  summarize(
    mean = mean(Chronological_Age),
    sd = sd(Chronological_Age, na.rm = TRUE),
    min = min(Chronological_Age),
    max = max(Chronological_Age)
  )

ss_primary %>%
  dplyr::summarize(
    n = sum(!is.na(Horvath_MethylClock)),
    median = median(Horvath_MethylClock, na.rm = TRUE),
    q1 = quantile(Horvath_MethylClock, 0.25, na.rm = TRUE),
    q3 = quantile(Horvath_MethylClock, 0.75, na.rm = TRUE),
    min = min(Horvath_MethylClock, na.rm = TRUE),
    max = max(Horvath_MethylClock, na.rm = TRUE),
    .groups = "drop"
  )
ss_primary %>%
  dplyr::summarize(
    n = sum(!is.na(Chronological_Age)),
    median = median(Chronological_Age, na.rm = TRUE),
    q1 = quantile(Chronological_Age, 0.25, na.rm = TRUE),
    q3 = quantile(Chronological_Age, 0.75, na.rm = TRUE),
    min = min(Chronological_Age, na.rm = TRUE),
    max = max(Chronological_Age, na.rm = TRUE),
    .groups = "drop"
  )

ss %>% filter(Primary_Include_In_Analysis == 1, Batch == "REF") %>% select(Histology_General) %>% table()
ss %>% filter(Lymph_Node == "T", Batch == "REF") %>% select(Histology_General) %>% table()
ss %>% filter(Primary_Include_In_Analysis == 1, Batch == "VAL") %>% select(Histology_General) %>% table()
ss %>% filter(Lymph_Node == "T", Batch == "VAL") %>% select(Histology_General) %>% table()

ss %>% filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>% select(Clinical_Invasiveness) %>% table()
ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>% select(Clinical_Invasiveness) %>% table()
ss %>% filter(Batch == "REF", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(Clinical_Invasiveness) %>% table()
ss %>% filter(Batch == "VAL", LN_Include_In_Analysis == 1, , Lymph_Node == "T") %>% select(Clinical_Invasiveness) %>% table()


ss %>% filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>% select(Tissue_Type) %>% table()
ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>% select(Tissue_Type) %>% table()
ss %>% filter(Batch == "REF", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(Tissue_Type) %>% table()
ss %>% filter(Batch == "VAL", LN_Include_In_Analysis == 1, , Lymph_Node == "T") %>% select(Tissue_Type) %>% table()


ss %>% filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>% select(Clinical_Confidence) %>% table()
ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>% select(Clinical_Confidence) %>% table()
ss %>% filter(Batch == "REF", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(Clinical_Confidence) %>% table()
ss %>% filter(Batch == "VAL", LN_Include_In_Analysis == 1, , Lymph_Node == "T") %>% select(Clinical_Confidence) %>% table()

ss %>% filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>% select(Driver_Group) %>% table()
ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>% select(Driver_Group) %>% table()
ss %>% filter(Batch == "REF", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(Driver_Group) %>% table()
ss %>% filter(Batch == "VAL", LN_Include_In_Analysis == 1, , Lymph_Node == "T") %>% select(Driver_Group) %>% table()

ss %>% filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>% select(T) %>% table()
ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>% select(T) %>% table()
ss %>% filter(Batch == "REF", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(T) %>% table()
ss %>% filter(Batch == "VAL", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(T) %>% table()

ss %>% filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>% select(N) %>% table()
ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>% select(N) %>% table()
ss %>% filter(Batch == "REF", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(N) %>% table()
ss %>% filter(Batch == "VAL", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(N) %>% table()

ss %>% filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>% select(M) %>% table()
ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>% select(M) %>% table()
ss %>% filter(Batch == "REF", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(M) %>% table()
ss %>% filter(Batch == "VAL", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(M) %>% table()

ss %>% filter(Batch == "REF", Primary_Include_In_Analysis == 1, M == "M1") %>% select(N) %>% table()
ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == 1,  M == "M1") %>% select(N) %>% table()
ss %>% filter(Batch == "REF", LN_Include_In_Analysis == 1, Lymph_Node == "T") %>% select(N) %>% table()
ss %>% filter(Batch == "VAL", LN_Include_In_Analysis == 1, , Lymph_Node == "T") %>% select(N) %>% table()



tab_ref <- ss %>%
  dplyr::filter(Batch == "REF", Primary_Include_In_Analysis == 1) %>%
  dplyr::pull(Driver_Group) %>%
  table()
cbind(
  Count = tab_ref,
  Proportion = prop.table(tab_ref)
)

tab_ref <- ss %>%
  dplyr::filter(Batch == "VAL", Primary_Include_In_Analysis == 1) %>%
  dplyr::pull(Driver_Group) %>%
  table()
cbind(
  Count = tab_ref,
  Proportion = prop.table(tab_ref)
)

# PSR 
ss %>%
  dplyr::filter(Batch == "REF", LN_Include_In_Analysis == "1") %>%
  dplyr::summarize(
    mean = mean(Probe_Success_Rate, na.rm = TRUE),
    sd   = sd(Probe_Success_Rate, na.rm = TRUE),
    min  = min(Probe_Success_Rate, na.rm = TRUE),
    q1   = quantile(Probe_Success_Rate, 0.25, na.rm = TRUE),
    q3   = quantile(Probe_Success_Rate, 0.75, na.rm = TRUE),
    max  = max(Probe_Success_Rate, na.rm = TRUE)
  )
ss %>%
  dplyr::filter(Batch == "VAL", LN_Include_In_Analysis == "1") %>%
  dplyr::summarize(
    mean = mean(Probe_Success_Rate, na.rm = TRUE),
    sd   = sd(Probe_Success_Rate, na.rm = TRUE),
    min  = min(Probe_Success_Rate, na.rm = TRUE),
    q1   = quantile(Probe_Success_Rate, 0.25, na.rm = TRUE),
    q3   = quantile(Probe_Success_Rate, 0.75, na.rm = TRUE),
    max  = max(Probe_Success_Rate, na.rm = TRUE)
  )


# ACCURACIES FROM FINAL MODEL for validation cohort
inv <- ss %>% filter(Batch == "VAL", Primary_Include_In_Analysis == "1") #85
dri <- inv %>% filter(Driver_Group != "Indeterminate") #78
dim(inv)
dim(dri)

# invasiveness wrong
wrong_inv <- inv %>% filter(Fmodel_Invasiveness_Accuracy == 1)
wrong_inv %>% select(Tissue_Type) %>% table()
wrong_inv %>% select(Clinical_Invasiveness) %>% table()
wrong_inv %>% filter(Clinical_Invasiveness == "High") %>%
  select(Driver_Group) %>% table()
wrong_inv %>% filter(Clinical_Invasiveness == "Low") %>%
  select(Driver_Group) %>% table()

right_inv <- inv %>% filter(Fmodel_Invasiveness_Accuracy == 0)
right_inv %>% select(Tissue_Type) %>% table()

wrong_inv %>% select(Clinical_Invasiveness) %>% table()
wrong_inv %>%
  summarize(
    accuracy = mean(IC_EpiDISH, na.rm = TRUE),
    sd = sd(as.numeric(IC_EpiDISH), na.rm = TRUE),
    min = min(IC_EpiDISH, na.rm = TRUE),
    max = max(IC_EpiDISH, na.rm = TRUE)
  )
right_inv %>%
  summarize(
    accuracy = mean(IC_EpiDISH, na.rm = TRUE),
    sd = sd(as.numeric(IC_EpiDISH), na.rm = TRUE),
    min = min(IC_EpiDISH, na.rm = TRUE),
    max = max(IC_EpiDISH, na.rm = TRUE)
  )

# driver wrong
dri %>% select(Fmodel_Driver_Accuracy) %>% table()
wrong_dri <- dri %>% filter(Fmodel_Driver_Accuracy == 1)
wrong_dri %>% select(Tissue_Type) %>% table()

right_dri <- dri %>% filter(Fmodel_Driver_Accuracy == 0)
right_dri %>% select(Tissue_Type) %>% table()

wrong_dri %>% select(Driver_Group) %>% table()
wrong_dri %>%
  summarize(
    accuracy = mean(IC_EpiDISH, na.rm = TRUE),
    sd = sd(as.numeric(IC_EpiDISH), na.rm = TRUE),
    min = min(IC_EpiDISH, na.rm = TRUE),
    max = max(IC_EpiDISH, na.rm = TRUE)
  )
right_dri %>%
  summarize(
    accuracy = mean(IC_EpiDISH, na.rm = TRUE),
    sd = sd(as.numeric(IC_EpiDISH), na.rm = TRUE),
    min = min(IC_EpiDISH, na.rm = TRUE),
    max = max(IC_EpiDISH, na.rm = TRUE)
  )

# ============================================================
# ADULT STATS
# ============================================================
ss <- as.data.frame(read_excel(file.path("ss", ADT_META)))

length(unique(ss$Sample_ID)) # 499
length(unique(ss$Patient_ID)) # 449 

ss %>% filter(Lymph_Node == "F") %>%
  summarize(
    mean = mean(Chronological_Age),
    sd = sd(Chronological_Age, na.rm = TRUE),
    min = min(Chronological_Age),
    max = max(Chronological_Age)
  )

ss %>% filter(Lymph_Node == "F") %>% select(Histology_General) %>% table()
ss %>% filter(Lymph_Node == "F") %>% select(Clinical_Invasiveness) %>% table()
ss %>% filter(Lymph_Node == "F") %>% select(Driver_Group) %>% table()

tab_ref <- ss %>%
  dplyr::filter(Lymph_Node == "F") %>%
  dplyr::pull(Driver_Group) %>%
  table()
cbind(
  Count = tab_ref,
  Proportion = prop.table(tab_ref)
)

ss %>% filter(Lymph_Node == "F") %>% 
  dplyr::summarize(
    n = sum(!is.na(Chronological_Age)),
    median = median(Chronological_Age, na.rm = TRUE),
    q1 = quantile(Chronological_Age, 0.25, na.rm = TRUE),
    q3 = quantile(Chronological_Age, 0.75, na.rm = TRUE),
    min = min(Chronological_Age, na.rm = TRUE),
    max = max(Chronological_Age, na.rm = TRUE),
    .groups = "drop"
  )

ss %>% filter(Lymph_Node == "F") %>% 
  dplyr::summarize(
    n = sum(!is.na(Horvath_MethylClock)),
    median = median(Horvath_MethylClock, na.rm = TRUE),
    q1 = quantile(Horvath_MethylClock, 0.25, na.rm = TRUE),
    q3 = quantile(Horvath_MethylClock, 0.75, na.rm = TRUE),
    min = min(Horvath_MethylClock, na.rm = TRUE),
    max = max(Horvath_MethylClock, na.rm = TRUE),
    .groups = "drop"
  )


