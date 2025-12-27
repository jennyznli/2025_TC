# ============================================================
# LOOCV INVASIVENESS ANALYSIS
# ============================================================
source("config.R")

CL_DIR <- file.path("data", "loocv_invasiveness_model")
MODEL_DIR <- file.path(CL_DIR, "models")
FEAT_DIR <- file.path(CL_DIR, "features")
SUM_DIR <- file.path(CL_DIR, "summary")

# ========================
# LOAD IN DATA
# ========================
pred_results <- read.csv(file.path(SUM_DIR, "loocv_all_predictions.csv"))
conf_matrix <- readRDS(file.path(SUM_DIR, "loocv_conf_matrix.rds"))
overall_metrics <- readRDS(file.path(SUM_DIR, "loocv_overall_metrics.rds"))
model_list <- readRDS(file.path(SUM_DIR, "loocv_model_list.rds"))
roc_obj <- readRDS(file.path(SUM_DIR, "loocv_roc_curve.rds"))
imp_results <- readRDS(file.path(SUM_DIR, "loocv_imp_results.rds"))

print(overall_metrics)
# Accuracy Balanced_Accuracy     Kappa       AUC Precision    Recall Specificity       F1
# Balanced Accuracy 0.8390805         0.8139205 0.6446908 0.8633523 0.8474576 0.9090909     0.71875 0.877193

print(conf_matrix)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction High Low
# High   50   9
# Low     5  23
# 
# Accuracy : 0.8391          
# 95% CI : (0.7448, 0.9091)
# No Information Rate : 0.6322          
# P-Value [Acc > NIR] : 1.892e-05       
# 
# Kappa : 0.6447          
# 
# Mcnemar's Test P-Value : 0.4227          
#                                           
#             Sensitivity : 0.9091          
#             Specificity : 0.7188          
#          Pos Pred Value : 0.8475          
#          Neg Pred Value : 0.8214          
#              Prevalence : 0.6322          
#          Detection Rate : 0.5747          
#    Detection Prevalence : 0.6782          
#       Balanced Accuracy : 0.8139          
#                                           
#        'Positive' Class : High      
       

# ========================
# JOIN iNTO SS 
# ========================
# check if agree w/ original LOOCV 
# ss <- read_excel(file.path("ss", PED_META)) %>% 
#   filter(Primary_Include_In_Analysis == 1, Reference == 1)
# ss1 <- left_join(ss, pred_results, by="IDAT")
# ss1 <- ss1 %>%
#   dplyr::mutate(
#     LOOCV_Invasiveness_Predicted2 = Predicted,
#     LOOCV_Invasiveness_Accuracy2 = ifelse(Predicted == Clinical_Invasiveness, 0, 1),
#   ) %>% 
#   dplyr::select(-True_Label, -Predicted, -Accuracy, -Fold, -Prob_High, -Prob_Low)
# 
# ss1$LOOCV_Agree_Invasiveness <- ss1$LOOCV_Invasiveness_Predicted2 == ss1$LOOCV_Predicted_Invasiveness
# table(ss1$LOOCV_Invasiveness_Predicted2, ss1$Clinical_Invasiveness)
# write.csv(ss1, file.path("data", "temp.csv"), na = "")

### JOIN W SS
ss <- read_excel(file.path("ss", PED_META))
ss1 <- left_join(ss, pred_results, by="IDAT")
ss1 <- ss1 %>%
  dplyr::mutate(
    LOOCV_Invasiveness_Prediction  = Predicted,
    LOOCV_Invasiveness_Accuracy   = ifelse(Predicted == Clinical_Invasiveness, 0, 1),
    LOOCV_Invasiveness_Confidence = pmax(Prob_High, Prob_Low)
  ) %>% 
  dplyr::select(-True_Label, -Predicted, -Accuracy, -Fold, -Prob_High, -Prob_Low)
write_xlsx(ss1, file.path("ss", "pediatric_thyroid_master.xlsx"))

# ========================
# ROC CURVE
# ========================
roc_data <- data.frame(
  Specificity = 1 - roc_obj$specificities,
  Sensitivity = roc_obj$sensitivities
)

p <- ggplot(roc_data, aes(x = Specificity, y = Sensitivity)) +
  geom_line(color = "#f8766dff", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  annotate("text", x = 0.75, y = 0.25,
           label = paste("AUC =", round(as.numeric(overall_metrics$AUC), 3)),
           size = 5, color = "#f8766dff") +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_equal()

pdf(file.path("figures", "loocv_invasiveness_roc.pdf"), height = 4, width = 4)
print(p)
dev.off()



make_plot <- function(df, acc_col, conf_col, title = NULL) {
  stopifnot(acc_col %in% colnames(df))
  stopifnot(conf_col %in% colnames(df))
  
  acc_raw <- df[[acc_col]]
  conf_raw <- df[[conf_col]]
  
  # ---- TYPE CHECKS ----
  if (!is.numeric(conf_raw)) {
    suppressWarnings(conf_num <- as.numeric(conf_raw))
    if (all(is.na(conf_num))) {
      stop(conf_col, " could not be coerced to numeric.")
    }
  } else {
    conf_num <- conf_raw
  }
  
  # Accuracy column validation
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
  
  # ---- FINAL SANITY CHECK ----
  if (!all(c("Correct","Incorrect") %in% df2$acc_bin)) {
    stop("Both Correct and Incorrect groups must be present for ", title)
  }
  
  # ---- ONE-SIDED WILCOXON TEST ----
  # Testing if Incorrect predictions have LOWER confidence than Correct
  wilcox_res <- wilcox.test(
    conf_val ~ acc_bin, 
    data = df2,
    alternative = "greater",  # Tests if Correct > Incorrect
    exact = FALSE
  )
  
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
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

ss_ref <- read_excel(file.path("ss", PED_META)) %>%
  filter(Include_In_Analysis == 1, Batch == "REF", Lymph_Node == "F")

p1 <- make_plot(
  df = ss_ref,
  acc_col = "LOOCV_Invasiveness_Accuracy",
  conf_col = "LOOCV_Invasiveness_Confidence",
  title = "Reference – Invasiveness"
)

ss_val <- read_excel(file.path("ss", PED_META)) %>%
  filter(Include_In_Analysis == 1, Batch == "VAL", Lymph_Node == "F")

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

ggsave(
  filename = file.path("figures", "all_accuracy_vs_confidence.pdf"),
  plot = final_plot,
  width = 10,
  height = 3
)


