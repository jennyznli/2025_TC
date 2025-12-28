# ============================================================
#   This script contains the analysis for final joint
#   invasiveness & driver classifier, including
#   joint probe importance scatter plot.
# ============================================================
BATCH_SIZE <- 10000
N_TREES <- 1000
N_FEATURES <- 3000

source("config.R")

R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
DRIV_DIR <- file.path(dir, "final_driver_model")
INV_DIR <- file.path(dir, "final_invasiveness_model")

# ========================
# SCATTER PLOT
# ========================
head <- 30000

imp_results_driver <- readRDS(file.path(DRIV_DIR, "fold_importance.rds"))
imp_results_invasive <- readRDS(file.path(INV_DIR, "fold_importance.rds"))

# Get top features
top_driver_features <- imp_results_driver %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(head) %>%
  pull(Feature)
top_invasive_features <- imp_results_invasive %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(head) %>%
  pull(Feature)

# Find intersection and unique features
common_features <- intersect(top_driver_features, top_invasive_features)
driver_unique <- setdiff(top_driver_features, top_invasive_features)
invasive_unique <- setdiff(top_invasive_features, top_driver_features)

# Create a combined importance data frame
combined_importance <- bind_rows(
  imp_results_driver %>%
    filter(Feature %in% c(common_features, driver_unique, invasive_unique)) %>%
    group_by(Feature) %>%
    summarize(Mean_Importance = mean(MeanDecreaseAccuracy, na.rm = TRUE)) %>%
    mutate(Model = "Driver"),
  imp_results_invasive %>%
    filter(Feature %in% c(common_features, driver_unique, invasive_unique)) %>%
    group_by(Feature) %>%
    summarize(Mean_Importance = mean(MeanDecreaseAccuracy, na.rm = TRUE)) %>%
    mutate(Model = "Invasiveness")
)

# Tag features as common or unique
combined_importance <- combined_importance %>%
  mutate(Feature_Type = case_when(
    Feature %in% common_features ~ "Common",
    Feature %in% driver_unique ~ "Driver-specific",
    Feature %in% invasive_unique ~ "Invasiveness-specific"
  ))

# SCATTER PLOT
importance_wide <- combined_importance %>%
  pivot_wider(id_cols = Feature,
              names_from = Model,
              values_from = Mean_Importance) %>%
  mutate(Feature_Type = case_when(
    Feature %in% common_features ~ "Common",
    Feature %in% driver_unique ~ "Driver-specific",
    Feature %in% invasive_unique ~ "Invasiveness-specific"
  ))

# Replace NA with zeros for visualization
importance_wide <- importance_wide %>%
  mutate(
    Driver = replace_na(Driver, 0),
    Invasiveness = replace_na(Invasiveness, 0)
  )

p <- ggplot(importance_wide,
            aes(x = Driver, y = Invasiveness, color = Feature_Type)) +
  geom_point(alpha = 0.7, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = feature_colors, name = "Feature Type") +
  # Add labels for top features
  geom_text_repel(
    data = importance_wide %>%
      filter(Driver > quantile(Driver, 0.9) |
               Invasiveness > quantile(Invasiveness, 0.9)),
    aes(label = Feature),
    size = 3,
    max.overlaps = 15,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "gray50"
  ) +
  labs(
    x = "Importance for Driver Classification",
    y = "Importance for Invasiveness Classification"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "lightgray", size = 0.2),
    panel.grid.minor = element_line(color = "lightgray", size = 0.1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

pdf(file.path(fig_dir, "fmodel_joint_importance_scatter.pdf"),
    height = 5, width = 5)
print(p)
dev.off()
