# ============================================================
#   This script trains the final RF classifier for
#   predicting driver group in thyroid cancer samples
#   using DNA methylation data and tests on validation cohort. 
# ============================================================
BATCH_SIZE <- 10000
N_TREES <- 1000
N_FEATURES <- 3000

source("config.R")

model_suffix <- "driver"
base_name <- tools::file_path_sans_ext(PED_BETAS)

MODEL_DIR <- file.path(BASE_DIR, "data", paste0("fmodel_", model_suffix, "_", base_name, "_all_samples"))
TRAIN_DIR <- file.path(MODEL_DIR, "train")
TEST_DIR <- file.path(MODEL_DIR, "test")
dir.create(TRAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TEST_DIR, recursive = TRUE, showWarnings = FALSE)

# ========================
# LOAD DATA
# ========================
ss <- readxl::read_excel(file.path("ss", PED_META))
ss <- ss[!is.na(ss[[TARGET]]), ]
ss <- ss[ss[[TARGET]] != "Indeterminate", ]

betas <- readRDS(file.path("data", PED_BETAS))

# split into train and test 
ss_train <- ss %>%
  dplyr::filter(Batch == 'REF') %>% 
  filter(Lymph_Node == "F")
ss_test <- ss %>%
  filter(!(Batch == 'REF' & Lymph_Node == "F"))

samples_train <- intersect(colnames(betas), ss_train$IDAT)
betas_train <- betas[, samples_train]
ss_train <- ss_train[match(samples_train, ss_train$IDAT), ]
stopifnot(identical(colnames(betas_train), ss_train$IDAT))

samples_test <- intersect(colnames(betas), ss_test$IDAT)
betas_test <- betas[, samples_test]
ss_test <- ss_test[match(samples_test, ss_test$IDAT), ]
stopifnot(identical(colnames(betas_test), ss_test$IDAT))

train_labels <- factor(ss_train[[TARGET]])
test_labels  <- factor(ss_test[[TARGET]], levels = levels(train_labels))

# ========================
# FEATURE SELECTION
# ========================
n_probes <- nrow(betas_train)
n_batches <- ceiling(n_probes / BATCH_SIZE)
fold_importance <- data.frame()

for (s in 0:(n_batches - 1)) {
  
  start_idx <- s * BATCH_SIZE + 1
  end_idx <- min((s + 1) * BATCH_SIZE, n_probes)
  sbetas <- betas_train[start_idx:end_idx, ]
  
  model <- randomForest::randomForest(
    x = t(sbetas), 
    y = train_labels, 
    ntree = N_TREES, 
    importance = TRUE
  )
  
  importance <- data.frame(
    Feature = rownames(model$importance),
    MeanDecreaseAccuracy = model$importance[, "MeanDecreaseAccuracy"]
  )
  fold_importance <- rbind(fold_importance, importance)
}

ranked_features <- fold_importance %>% 
  dplyr::arrange(desc(MeanDecreaseAccuracy)) %>%
  dplyr::pull(Feature)

saveRDS(ranked_features, file.path(TRAIN_DIR, "ranked_features.rds"))
saveRDS(fold_importance, file.path(TRAIN_DIR, "fold_importance.rds"))

# ========================
# TRAIN FINAL MODEL
# ========================
sel_probes <- ranked_features[1:N_FEATURES]
train_betas_sel <- betas_train[sel_probes, ]

model <- randomForest::randomForest(
  x = t(train_betas_sel),
  y = train_labels,
  ntree = N_TREES,
  importance = TRUE
)

saveRDS(model, file.path(TRAIN_DIR, "final_model.rds"))

# ========================
# TEST PREDICTIONS
# ========================
test_betas_sel <- betas_test[sel_probes, ]

pred_test_class <- predict(model, t(test_betas_sel))
pred_test_prob <- predict(model, t(test_betas_sel), type = "prob")

predictions <- data.frame(
  IDAT = colnames(test_betas_sel),
  Predicted_Class = pred_test_class,
  pred_test_prob,  
  Max_Probability = apply(pred_test_prob, 1, max)
)

write.csv(predictions, file.path(TEST_DIR, "fmodel_pred_probs.csv"), row.names = FALSE)

# ========================
# SAVE RESULTS
# ========================
ss_full <- as.data.frame(readxl::read_excel(file.path("ss", PED_META)))
results <- dplyr::left_join(ss_full, predictions, by = "IDAT")
results$Correct_Prediction <- (results$Driver_Group == results$Predicted_Class)

writexl::write_xlsx(results, file.path("ss", "ss_fmodel_driver_predictions.xlsx"))

# ========================
# JOIN SPREADSHEET & METRICS
# ========================
results <- read.csv(file.path("data", "fmodel_driver_ped_betas_QCDPB_prc_all_samples", "test", "fmodel_pred_probs.csv"))

ss <- readxl::read_excel(file.path("ss", PED_META)) %>% filter(Driver_Group != "Indeterminate")
ss$Driver_Group <- as.factor(ss$Driver_Group)

ss <- left_join(ss, results, by="IDAT") %>% 
  mutate(Fmodel_Driver_Prediction = Predicted_Class,
         Fmodel_Driver_Confidence = Max_Probability, 
         Fmodel_Driver_Accuracy = ifelse(Fmodel_Driver_Prediction == Driver_Group, 0, 1)
  ) %>% 
  select(-DICER1, -Kinase.Fusion, -Ras.like, -BRAF.V600E, -Max_Probability, -Predicted_Class)

ss$Fmodel_Driver_Prediction <- factor(
  ss$Fmodel_Driver_Prediction,
  levels = levels(ss$Driver_Group)
)

write_xlsx(ss, file.path("ss", "ped_fmodel_driver_ss.xlsx"))

ss_val <- ss %>% 
  filter(Batch == "VAL", Lymph_Node == "F")

ss_ln <- ss %>% 
  filter(Lymph_Node == "T")

cm_val <- confusionMatrix(data = ss_val$Fmodel_Driver_Prediction, reference = ss_val$Driver_Group)
print(cm_val)

# Confusion Matrix and Statistics
# 
# Reference
# Prediction      BRAF V600E DICER1 Kinase Fusion Ras-like
# BRAF V600E            14      0             4        0
# DICER1                 0      1             0        0
# Kinase Fusion          5      0            39        0
# Ras-like               1      1             2        9
# 
# Overall Statistics
# 
# Accuracy : 0.8289          
# 95% CI : (0.7253, 0.9057)
# No Information Rate : 0.5921          
# P-Value [Acc > NIR] : 8.466e-06       
# 
# Kappa : 0.7021          
# 
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: BRAF V600E Class: DICER1 Class: Kinase Fusion Class: Ras-like
# Sensitivity                     0.7000       0.50000               0.8667          1.0000
# Specificity                     0.9286       1.00000               0.8387          0.9403
# Pos Pred Value                  0.7778       1.00000               0.8864          0.6923
# Neg Pred Value                  0.8966       0.98667               0.8125          1.0000
# Prevalence                      0.2632       0.02632               0.5921          0.1184
# Detection Rate                  0.1842       0.01316               0.5132          0.1184
# Detection Prevalence            0.2368       0.01316               0.5789          0.1711
# Balanced Accuracy               0.8143       0.75000               0.8527          0.9701


wrong <- ss_val %>% filter(Fmodel_Driver_Accuracy == 1, Lymph_Node == "F", Driver_Group != "Indeterminate")

### LYMPH NODE EVALUATION 
cm_ln <- confusionMatrix(data = ss_ln$Fmodel_Driver_Prediction, reference = ss_ln$Driver_Group)
print(cm_ln)
# Reference
# Prediction      BRAF V600E DICER1 Kinase Fusion Ras-like
# BRAF V600E             2      0             0        0
# DICER1                 0      0             0        0
# Kinase Fusion          0      0            10        0
# Ras-like               0      1             0        1

# ========================
# ONE V ALL ROC CURVES
# ========================
ss_val <- ss_val %>% left_join(results, by="IDAT")

ss_val <- ss_val %>%
  rename(
    "BRAF.V600E" = "BRAF V600E",
    "Ras.like" = "Ras-like",
    "Kinase.Fusion" = "Kinase Fusion"
  )

driver_classes <- levels(ss_val$Driver_Group)
roc_list_val <- list()
roc_df_val <- data.frame()

for (cls in driver_classes) {
  roc_obj <- roc(response  = ss_val$Driver_Group == cls, predictor = ss_val[[cls]], quiet = TRUE)
  roc_list_val[[cls]] <- roc_obj
  
  roc_df_val <- rbind(
    roc_df_val,
    data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Class = cls,
      AUC = as.numeric(roc_obj$auc)
    )
  )
}

auc_labels <- roc_df_val %>%
  group_by(Class) %>%
  summarize(AUC = unique(AUC)) %>%
  mutate(label = paste0(Class, " (AUC = ", round(AUC, 3), ")"))

p <- ggplot(roc_df_val, aes(FPR, TPR, color = Class)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = driver_colors, labels = auc_labels$label) +
  coord_equal() +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "False Positive Rate", y = "True Positive Rate")

pdf(file.path("figures", "val_driver_roc.pdf"), width = 5, height = 5)
print(p)
dev.off()


# ========================
# ENRICHMENT DOT PLOT - MOST IMPORTANT
# ========================
imp <- readRDS(file.path("data/fmodel_driver_ped_betas_QCDPB_prc_all_samples/train/ranked_features.rds"))

sel_probes3k <- imp %>% head(3000)

x <- testEnrichment(sel_probes3k, "TFBSconsensus", platform = "EPIC")
pdf(file.path("figures", "fmodel_driver_enrichment_tfbs.pdf"),width = 2.8, height = 3.9, onefile=FALSE)
plotDotFDR(x, n_min = 30, n_max = 30)
dev.off()


# ========================
# ONE V ALL SHAP MODEL S
# ========================
train_betas_sel <- readRDS(file.path("data", "train_betas_sel.rds"))
unique_classes <- levels(train_labels)

for (class_name in unique_classes) {
  binary_target <- as.numeric(train_labels == class_name)
  train_data <- t(train_betas_sel)
  
  binary_model <- randomForest(
    x = train_data,
    y = as.factor(binary_target),
    mtry = 2,
    nodesize = 5,
    importance = TRUE
  )
  clean <- clean_label(class_name)
  saveRDS(binary_model, file.path(MODEL_DIR, paste0(clean, "_shap_model.rds")))
  
  unified_model <- randomForest.unify(binary_model, train_data)
  shap <- treeshap(unified_model, train_data, interactions = TRUE)
  shp <- shapviz::shapviz(shap, X = train_data)
  
  saveRDS(shap, file.path(SHAP_DIR, paste0(clean, "_shap.rds")))
  saveRDS(shp, file.path(SHAP_DIR, paste0(clean, "_shpviz.rds")))
  
  # most important features
  most_imp <- sort(colMeans(abs(shap$shaps)), decreasing = TRUE)
  # saveRDS(most_imp, file.path(SHAP_DIR, paste0(clean, "_shap_importance.rds")))
  write.csv(data.frame(Feature = names(most_imp),
                       SHAP_Value = most_imp),
            file.path(SHAP_DIR, paste0(clean, "_shap_importance.csv")))
}

# ========================
# SHAP PLOTS
# ========================
for (class_name in unique_classes) {
  clean <- clean_label(class_name)
  
  # read in
  shp <- readRDS(file.path(SHAP_DIR, paste0(clean, "_shpviz.rds")))
  shap <- readRDS(file.path(SHAP_DIR, paste0(clean, "_shap.rds")))
  
  # beeswarm
  p <- sv_importance.shapviz(shp, kind = "beeswarm", max_display = 30)
  pdf(file.path(FIG_DIR, paste0(clean, "_ped85_beeswarm.pdf")), height = 6, width = 5)
  plot(p)
  dev.off()
  
  # contribution plot
  # p <- plot_contribution(shap, obs = 1, max_vars = 30)
  # pdf(file.path(fig_dir, paste0(clean, "_ped85_contribution.pdf")), height = 7, width = 10)
  # plot(p)
  # dev.off()
}
