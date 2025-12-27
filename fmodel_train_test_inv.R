# ============================================================
#   FIG 4.
#   This script trains the final RF classifier for
#   predicting clinical invasiveness in thyroid cancer samples
#   using DNA methylation data.
# ============================================================
TARGET <- "Clinical_Invasiveness"
BATCH_SIZE <- 10000
N_TREES <- 1000
N_FEATURES <- 3000

source("config.R")

model_suffix <- if (TARGET == "Clinical_Invasiveness") "invasiveness" else "driver"
base_name <- tools::file_path_sans_ext(PED_BETAS)

MODEL_DIR <- file.path("models", paste0("fmodel_", model_suffix, "_", base_name, "_all_samples"))
TRAIN_DIR <- file.path(MODEL_DIR, "train")
TEST_DIR <- file.path(MODEL_DIR, "test")
dir.create(TRAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TEST_DIR, recursive = TRUE, showWarnings = FALSE)

# ========================
# LOAD DATA
# ========================
ss <- readxl::read_excel(file.path("ss", PED_META))
ss <- ss[!is.na(ss[[TARGET]]), ]

betas <- readRDS(file.path("data", PED_BETAS))

# Split into train (Primary batch) and test (all others)
ss_train <- ss %>%
  dplyr::filter(Batch == 'REF') %>% 
  filter(Primary_Include_In_Analysis == "1")

ss_test <- ss %>%
  filter(!(Batch == 'REF' & Primary_Include_In_Analysis == "1"))

# Match samples
samples_train <- intersect(colnames(betas), ss_train$IDAT)
betas_train <- betas[, samples_train]
ss_train <- ss_train[match(samples_train, ss_train$IDAT), ]
stopifnot(identical(colnames(betas_train), ss_train$IDAT))

samples_test <- intersect(colnames(betas), ss_test$IDAT)
betas_test <- betas[, samples_test]
ss_test <- ss_test[match(samples_test, ss_test$IDAT), ]
stopifnot(identical(colnames(betas_test), ss_test$IDAT))

train_labels <- as.factor(ss_train[[TARGET]])
test_labels <- as.factor(ss_test[[TARGET]])

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

ss_full <- as.data.frame(readxl::read_excel(file.path("ss", PED_META)))
results <- dplyr::left_join(ss_full, predictions, by = "IDAT")
results$Accuracy <- (results$Clinical_Invasiveness == results$Predicted_Class)

writexl::write_xlsx(results, file.path("ss", "ss_fmodel_invasiveness_predictions.xlsx"))

# ========================
# JOIN TO SS
# ========================
results <- read.csv(file.path("data", "fmodel_invasiveness_ped_betas_QCDPB_prc_all_samples", "test", "fmodel_pred_probs.csv"))

ss_val <- readxl::read_excel(file.path("ss", PED_META)) %>% filter(Primary_Include_In_Analysis == "1", Batch == "VAL")
ss_val$Clinical_Invasiveness <- factor(ss_val$Clinical_Invasiveness)

ss_val <- left_join(ss_val, results, by="IDAT") %>% 
  mutate(Fmodel_Invasiveness_Prediction = Predicted_Class,
         Fmodel_Invasiveness_Confidence = Max_Probability, 
         Fmodel_Invasiveness_Accuracy = ifelse(Fmodel_Invasiveness_Prediction == Clinical_Invasiveness, 0, 1),
         Fmodel_Invasiveness_High_Prob = High, 
         Fmodel_Invasiveness_Low_Prob = Low
         ) %>% 
  select(-Low, -High, -Max_Probability, -Predicted_Class)

ss_val$Fmodel_Invasiveness_Prediction <- factor(
  ss_val$Fmodel_Invasiveness_Prediction,
  levels = levels(ss_val$Clinical_Invasiveness)
)

write_xlsx(ss_val, file.path("ss", "ped_fmodel_inv_ss.xlsx"))

conf_mat <- confusionMatrix(
  data = ss_val$Fmodel_Invasiveness_Prediction,
  reference = ss_val$Clinical_Invasiveness
)

print(conf_mat)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction High Low
# High   55  16
# Low     3   9
# 
# Accuracy : 0.7711          
# 95% CI : (0.6658, 0.8562)
# No Information Rate : 0.6988          
# P-Value [Acc > NIR] : 0.091803        
# 
# Kappa : 0.3618          
# 
# Mcnemar's Test P-Value : 0.005905        
#                                           
#             Sensitivity : 0.9483          
#             Specificity : 0.3600          
#          Pos Pred Value : 0.7746          
#          Neg Pred Value : 0.7500          
#              Prevalence : 0.6988          
#          Detection Rate : 0.6627          
#    Detection Prevalence : 0.8554          
#       Balanced Accuracy : 0.6541          
#                                           
#        'Positive' Class : High           

wrong <- ss_val %>% filter(Fmodel_Invasiveness_Accuracy == 1) %>% arrange(Clinical_Invasiveness) 
wrong %>% select(Clinical_Invasiveness) %>% table()

# ========================
# ROC CURVE (VALIDATION)
# ========================
roc_obj <- roc(
  response  = ss_val$Clinical_Invasiveness,
  predictor = ss_val$Fmodel_Invasiveness_High_Prob,
  levels    = c("Low", "High"),
  direction = "<"
)

auc_val <- as.numeric(auc(roc_obj))
print(auc_val)

roc_data <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#f8766dff", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "gray70") +
  annotate(
    "text",
    x = 0.65, y = 0.25,
    label = paste0("AUC = ", round(auc_val, 3)),
    size = 4.5,
    color = "#f8766dff"
  ) +
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10)
  )

pdf(file.path("figures", "val_invasiveness_roc.pdf"),
    height = 4, width = 4)
print(p)
dev.off()

# ========================
# ENRICHMENT PLOT
# ========================
imp <- readRDS(file.path("data/fmodel_invasiveness_ped_betas_QCDPB_prc_all_samples/train/ranked_features.rds"))
sel_probes3k <- imp %>% head(3000)

x <- testEnrichment(sel_probes3k, "TFBSconsensus", platform = "EPIC")
pdf(file.path("figures", "fmodel_invasiveness_enrichment_tfbs.pdf"),width = 2.8, height = 3.9, onefile=FALSE)
plotDotFDR(x, n_min = 30, n_max = 30)
dev.off()


