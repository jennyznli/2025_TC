# ============================================================
# Performs joint epigenetic age acceleration analysis. 
# ============================================================
source("config.R")

# ========================
# LOAD DATA
# ========================
betas_ped <- readRDS(file.path("data", "ped_betas_QCDPB_prc.rds"))
betas_adt <- readRDS(file.path("data", "adt_betas_QCDPB_pr.rds"))

ss_ped <- read_excel(file.path("ss", PED_META))
ss_adt <- read_excel(file.path("ss", ADT_META))

# ========================
# METHYL CLOCK INFERENCE
# ========================
cpgs.missing.ped <- checkClocks(betas_ped)
cpgs.missing.adt <- checkClocks(betas_adt)

age_ped <- DNAmAge(betas_ped)
age_adt <- DNAmAge(betas_adt)

ss_adt <- left_join(ss_adt, age_adt, by = "IDAT")
ss_ped <- left_join(ss_ped, age_adt, by = "IDAT")

write.csv(ss_ped, file.path("data", "ped_betas_QCDPB_pr_methylclock.csv"), row.names=FALSE)
write.csv(ss_adt, file.path("data", "adt_betas_QCDPB_pr_methylclock.csv"), row.names=FALSE)

# ========================
# EPIGENETIC AGE ACCELERATION
# ========================
ss <- read_xlsx(file.path("ss", JNT_META)) %>% filter(Reference == "1", Lymph_Node == "F")

ss$Chronological_Age <- as.numeric(ss$Chronological_Age)
ss$Epigenetic_Age <- as.numeric(ss$Epigenetic_Age)
ss <- ss %>%
  mutate(EAA = residuals(lm(Epigenetic_Age ~ Chronological_Age + Sex + IC_EpiDISH))) %>%
  filter(Clinical_Invasiveness %in% c("High", "Low"), !is.na(EAA))

# ========================
# CORRELATIONS
# ========================
invasiveness_levels <- unique(ss$Clinical_Invasiveness[!is.na(ss$Clinical_Invasiveness)])

cor_results <- list()
for (inv in invasiveness_levels) {
  subset_data <- ss %>% filter(Clinical_Invasiveness == inv)
  cor_inv <- cor.test(subset_data$Chronological_Age, subset_data$Epigenetic_Age,
                      method = "pearson", use = "complete.obs")
  cor_results[[inv]] <- cor_inv
  cat(sprintf("%s: r = %.3f, p-value = %.2e, n = %d\n",
              inv, cor_inv$estimate, cor_inv$p.value,
              sum(complete.cases(subset_data$Chronological_Age, subset_data$Epigenetic_Age))))
}
# Low: r = 0.726, p-value = 2.19e-43, n = 257
# High: r = 0.778, p-value = 9.27e-58, n = 279

# ========================
# SCATTER PLOT - CHRONOLOGICAL VS EPIGENETIC
# ========================
annotations <- data.frame(
  Clinical_Invasiveness = names(cor_results),
  r = sapply(cor_results, function(x) x$estimate),
  label = sapply(cor_results, function(x) sprintf("r = %.3f", x$estimate))
)

annotations$x <- c(min(ss$Chronological_Age, na.rm = TRUE) + 5, 
                   max(ss$Chronological_Age, na.rm = TRUE) - 15)
annotations$y <- c(max(ss$Epigenetic_Age, na.rm = TRUE) - 5,
                   max(ss$Epigenetic_Age, na.rm = TRUE) - 15)

p <- ggplot(ss, aes(x = Chronological_Age, y = Epigenetic_Age)) +
  geom_point(aes(color = Clinical_Invasiveness, shape = Source), size = 1) +
  scale_color_manual(values = invasiveness_colors) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", size = 0.3) +
  scale_shape_manual(values = source_shapes) +
  geom_text(data = annotations, 
            aes(x = x, y = y, label = label, color = Clinical_Invasiveness),
            hjust = 0, vjust = 1, size = 3.5, fontface = "bold", show.legend = FALSE) +
  labs(x = "Chronological Age", y = "Predicted Age") +
  theme_minimal()

ggsave(file.path("figures", "joint_actual_predicted_age_inv_with_corr.pdf"),
       plot = p, width = 5, height = 4)


# ========================
# BOX PLOT - INVASIVENESS
# ========================
format_pval <- function(p) {
  if(p < 0.001) return("p < 0.001")
  return(sprintf("p = %.3f", p))
}
stats_list <- lapply(unique(ss$Source), function(src) {
  data_subset <- ss[ss$Source == src, ]
  test <- wilcox.test(EAA ~ Clinical_Invasiveness, data = data_subset)
  data.frame(
    Source = src,
    p.value = test$p.value,
    label = format_pval(test$p.value),
    y.position = max(data_subset$EAA, na.rm = TRUE) + 2
  )
})
stats_df <- do.call(rbind, stats_list)

p <- ggplot(ss, aes(x = Source, y = EAA, fill = Clinical_Invasiveness)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.3, size = 0.6) +
  scale_fill_manual(values = invasiveness_colors) +
  theme_minimal() +
  labs(x = "Cohort",
       y = "Age Acceleration (Years)") +
  theme(axis.text.x = element_text(hjust = 1)) +
  geom_text(data = stats_df,
            aes(x = Source, y = y.position, label = label),
            inherit.aes = FALSE,
            size = 3)
ggsave(file.path("figures", "age_acceleration_invasiveness_cohort.pdf"),
       plot = p, width=5.5, height=4.5)

print("Wilcoxon test results by cohort:")
print(stats_df)

# Source    p.value     label y.position
# 1    PED 0.03900595 p = 0.039   25.69672
# 2    ADT 0.12908868 p = 0.129   69.18577