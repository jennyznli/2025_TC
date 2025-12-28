# ============================================================
# Plots figures for LOOCV of invasiveness classifier.  
# ============================================================
source("config.R")

CL_DIR <- file.path("data", "loocv_invasiveness_model")
MODEL_DIR <- file.path(CL_DIR, "models")
FEAT_DIR <- file.path(CL_DIR, "features")
SUM_DIR <- file.path(CL_DIR, "summary")

# ========================
# MODEL ANALYSIS
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
       
# JOIN WITH SS
ss <- read_excel(file.path("ss", PED_META))
ss1 <- left_join(ss, pred_results, by="IDAT")
ss1 <- ss1 %>%
  dplyr::mutate(
    LOOCV_Invasiveness_Prediction  = Predicted,
    LOOCV_Invasiveness_Accuracy   = ifelse(Predicted == Clinical_Invasiveness, 0, 1),
    LOOCV_Invasiveness_Confidence = pmax(Prob_High, Prob_Low)
  ) %>% 
  dplyr::select(-True_Label, -Predicted, -Accuracy, -Fold, -Prob_High, -Prob_Low)
write_xlsx(ss1, file.path("ss", PED_META))

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
