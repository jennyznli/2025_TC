# ============================================================
# Plots figures for LOOCV of driver classifier.  
# ============================================================
source("config.R")

CL_DIR <- file.path("data", "loocv_driver_model")
MODEL_DIR <- file.path(CL_DIR, "models")
FEAT_DIR <- file.path(CL_DIR, "features")
SUM_DIR <- file.path(CL_DIR, "summary")

# ========================
# MODEL ANALYSIS
# ========================
overall_metrics <- readRDS(file.path(SUM_DIR, "loocv_overall_metrics.rds"))
class_metrics <- readRDS(file.path(SUM_DIR, "loocv_per_class_metrics.rds"))
pred_results <- read.csv(file.path(SUM_DIR, "loocv_all_predictions.csv"))
conf_matrix <- readRDS(file.path(SUM_DIR, "loocv_conf_matrix.rds"))
model_list <- readRDS(file.path(SUM_DIR, "loocv_model_list.rds"))
roc_obj <- readRDS(file.path(SUM_DIR, "loocv_roc_curves.rds"))
imp_results <- readRDS(file.path(SUM_DIR, "loocv_imp_results.rds"))

print(overall_metrics)
# Accuracy Overall_Kappa AUC_BRAF V600E AUC_DICER1 AUC_Kinase Fusion AUC_Ras-like  Mean_AUC Macro_Sensitivity
# Kappa 0.952381     0.9323399      0.9926686          1         0.9714706    0.9757813 0.9849801         0.9614305
# Macro_Specificity Macro_Precision Macro_Recall  Macro_F1 Macro_BalancedAccuracy
# Kappa         0.9831552       0.9572511    0.9614305 0.9591497              0.9722929

print(class_metrics)
# Class Sensitivity Specificity Precision    Recall        F1 BalancedAccuracy
# Sensitivity     BRAF V600E   0.9545455   0.9838710 0.9545455 0.9545455 0.9545455        0.9692082
# Sensitivity1        DICER1   1.0000000   1.0000000 1.0000000 1.0000000 1.0000000        1.0000000
# Sensitivity2 Kinase Fusion   0.9411765   0.9800000 0.9696970 0.9411765 0.9552239        0.9605882
# Sensitivity3      Ras-like   0.9500000   0.9687500 0.9047619 0.9500000 0.9268293        0.9593750
# 1            Macro-Average   0.9614305   0.9831552 0.9572511 0.9614305 0.9591497        0.9722929

print(conf_matrix)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction      BRAF V600E DICER1 Kinase Fusion Ras-like
# BRAF V600E            21      0             1        0
# DICER1                 0      8             0        0
# Kinase Fusion          0      0            32        1
# Ras-like               1      0             1       19
# 
# Overall Statistics
# 
# Accuracy : 0.9524          
# 95% CI : (0.8825, 0.9869)
# No Information Rate : 0.4048          
# P-Value [Acc > NIR] : < 2.2e-16       
# 
# Kappa : 0.9323          
# 
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: BRAF V600E Class: DICER1 Class: Kinase Fusion Class: Ras-like
# Sensitivity                     0.9545       1.00000               0.9412          0.9500
# Specificity                     0.9839       1.00000               0.9800          0.9688
# Pos Pred Value                  0.9545       1.00000               0.9697          0.9048
# Neg Pred Value                  0.9839       1.00000               0.9608          0.9841
# Prevalence                      0.2619       0.09524               0.4048          0.2381
# Detection Rate                  0.2500       0.09524               0.3810          0.2262
# Detection Prevalence            0.2619       0.09524               0.3929          0.2500
# Balanced Accuracy               0.9692       1.00000               0.9606          0.9594

# JOIN W SS
ss <- read_excel(file.path("ss", PED_META))
ss1 <- left_join(ss, pred_results, by="IDAT")
ss1 <- ss1 %>%
  dplyr::mutate(
    LOOCV_Driver_Prediction = Predicted,
    LOOCV_Driver_Accuracy   = ifelse(Predicted == Driver_Group, 0, 1),
    LOOCV_Driver_Confidence = pmax(Prob_BRAF.V600E, Prob_DICER1, Prob_Kinase.Fusion, Prob_Ras.like)
  ) %>% 
  dplyr::select(-True_Label, -Predicted, -Accuracy, -Fold, -Prob_BRAF.V600E, -Prob_DICER1, -Prob_Kinase.Fusion, -Prob_Ras.like)

write_xlsx(ss1, file.path("ss", PED_META))

# ========================
# ROC CURVES
# ========================
roc_data <- data.frame()
class_names <- names(roc_obj)
n_curves <- length(roc_obj)

for (i in 1:length(roc_obj)) {
  current_roc <- roc_obj[[i]]
  current_auc <- as.numeric(current_roc$auc)
  print(current_auc)
  temp_data <- data.frame(
    FPR = 1 - current_roc$specificities,
    TPR = current_roc$sensitivities,
    Class = class_names[i],
    AUC = current_auc
  )
  print(dim(temp_data))
  roc_data <- rbind(roc_data, temp_data)
}
dim(roc_data)

p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Class)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = driver_colors,
                     labels = paste0(class_names, " (AUC = ",
                                     round(tapply(roc_data$AUC, roc_data$Class, unique), 3), ")")) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  coord_equal()

pdf(file.path("figures", "loocv_driver_roc.pdf"),
    height = 5,
    width = 5)
print(p)
dev.off()

