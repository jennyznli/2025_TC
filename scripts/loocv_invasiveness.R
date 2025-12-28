# ============================================================
#   This script performs leave-one-out cross validation (LOOCV)
#   to develop and validate a random forest classifier for
#   predicting clinical invasiveness in thyroid cancer samples
#   using DNA methylation data.
# ============================================================
N_WORKERS <- 20
set.seed(123)

source(("../config.R"))

CL_DIR <- file.path(BASE_DIR, "data", "loocv_invasiveness_model")
MODEL_DIR <- file.path(CL_DIR, "models")
FEAT_DIR <- file.path(CL_DIR, "features")
SUM_DIR <- file.path(CL_DIR, "summary")
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FEAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUM_DIR, recursive = TRUE, showWarnings = FALSE)

models_list <- list()
pred_results <- data.frame()
imp_results <- data.frame()

# ========================
# FUNCTIONS
# ========================
subset_feature_selection <- function(betas, batch_size, n_feat, train_idx, labels, n_trees = 500, fold = NULL) {
    n_probes <- nrow(betas)
    n_batches <- ceiling(n_probes / batch_size)

    fold_models <- list()
    fold_importance <- data.frame()

    for (s in 0:(n_batches - 1)) {
        # Calculate indices for current batch
        start_idx <- s * batch_size + 1
        end_idx <- min((s + 1) * batch_size, n_probes)
        sbetas <- betas[start_idx:end_idx, , drop = FALSE]
        print(paste("Batch dimensions:", dim(sbetas)[1], "x", dim(sbetas)[2]))

        model <- randomForest(
            x = t(sbetas[, train_idx]),
            y = as.factor(labels[train_idx]),
            ntree = n_trees,
            importance = TRUE
        )
        fold_models[[s + 1]] <- model

        importance <- data.frame(
            Fold = if (is.null(fold)) NA else fold,
            Model = s,
            Feature = rownames(model$importance),
            MeanDecreaseAccuracy = model$importance[, "MeanDecreaseAccuracy"]
        )
        fold_importance <- rbind(fold_importance, importance)
    }

    features <- fold_importance %>%
        arrange(desc(MeanDecreaseAccuracy)) %>%
        head(n_feat) %>%
        pull(Feature)

    return(features)
}

process_sample <- function(i, all_samples, betas, labels) {
    print(paste("Processing sample", i, "of", length(all_samples)))

    # Create LOOCV fold: one test sample, all others for training
    valid_idx <- all_samples[i]
    train_idx <- all_samples[-i]

    # Feature selection
    sel_probes <- subset_feature_selection(
        betas = betas,
        batch_size = 10000,
        n_feat = 3000,
        train_idx = train_idx,
        labels = labels,
        fold = i
    )

    model <- randomForest(
        x = t(betas[sel_probes, train_idx]),
        y = labels[train_idx],
        ntree = 500,
        importance = TRUE
    )

    importance <- data.frame(
        Feature = rownames(model$importance),
        MeanDecreaseAccuracy = model$importance[, "MeanDecreaseAccuracy"],
        MeanDecreaseGini = model$importance[, "MeanDecreaseGini"]
    )

    # Add class-specific importance
    for (level in levels(labels)) {
        if (level %in% colnames(model$importance)) {
            importance[, level] <- model$importance[, level]
        }
    }

    importance$Fold <- i

    # Predictions
    pred_probs <- predict(model, t(betas[sel_probes, valid_idx, drop = FALSE]), type = "prob")
    pred_class <- predict(model, t(betas[sel_probes, valid_idx, drop = FALSE]))

    preds <- data.frame(
        IDAT = ss$IDAT[valid_idx],
        True_Label = labels[valid_idx],
        Predicted = pred_class,
        Accuracy = (pred_class == labels[valid_idx]),
        Fold = i
    )

    # Add probability columns for each class
    for (j in 1:length(levels(labels))) {
        level_name <- levels(labels)[j]
        # Create a safe column name that matches R's automatic name conversion
        safe_level_name <- gsub(" ", ".", gsub("-", ".", level_name))
        col_name <- paste0("Prob_", safe_level_name)
        preds[, col_name] <- pred_probs[, j]
        # Store the mapping between original level and column name
        if (j == 1) {
            preds$LevelToColMapping <- paste(level_name, "->", col_name)
        }
    }

    preds$LevelToColMapping <- NULL

    saveRDS(model, file.path(MODEL_DIR, paste0("model_sample", i, ".rds")))
    saveRDS(importance, file.path(FEAT_DIR, paste0("importance_sample", i, ".rds")))

    return(list(
        sample = i,
        prediction = pred_class,
        prediction_probs = pred_probs,
        true_label = labels[valid_idx],
        preds = preds,
        importance = importance,
        selected_features = sel_probes,
        model = model,
        status = "complete"
    ))
}

# ========================
# MAIN EXECUTION
# ========================
ss <- read_excel(file.path("ss", PED_META))
ss <- ss %>% filter(Lymph_Node == "F")
betas <- readRDS(file.path("data", PED_BETAS))

valid_ids <- intersect(ss$IDAT, colnames(betas))
ss <- ss[match(valid_ids, ss$IDAT), ]
betas <- betas[, valid_ids]

labels <- as.factor(ss$Clinical_Invasiveness)

all_samples <- 1:length(labels)
print(paste("Total samples for LOOCV:", length(all_samples)))

print("Starting LOOCV")
results <- bplapply(all_samples, function(i) process_sample(i, all_samples, betas, labels), BPPARAM = MulticoreParam(N_WORKERS))
print("Completed all samples")

# ========================
# ANALYSIS
# ========================
for (i in 1:length(results)) {
    models_list[[i]] <- results[[i]]$model
    pred_results <- rbind(pred_results, results[[i]]$preds)
    imp_results <- rbind(imp_results, results[[i]]$importance)
}

conf_matrix <- caret::confusionMatrix(
    data = as.factor(pred_results$Predicted),
    reference = as.factor(pred_results$True_Label)
)

accuracy <- mean(pred_results$Accuracy)

positive_class <- "High"
prob_col_name <- paste0("Prob_", positive_class)

roc_obj <- pROC::roc(
    response = pred_results$True_Label,
    predictor = pred_results[[prob_col_name]],
    levels = levels(as.factor(pred_results$True_Label))
)

auc_value <- pROC::auc(roc_obj)

overall_metrics <- data.frame(
    Accuracy = accuracy,
    Balanced_Accuracy = conf_matrix$byClass["Balanced Accuracy"],
    Kappa = conf_matrix$overall["Kappa"],
    AUC = as.numeric(auc_value),
    Precision = conf_matrix$byClass["Pos Pred Value"],
    Recall = conf_matrix$byClass["Sensitivity"],
    Specificity = conf_matrix$byClass["Specificity"],
    F1 = conf_matrix$byClass["F1"]
)

# ========================
# SAVE FINAL RESULTS
# ========================
write.csv(pred_results, file.path(SUM_DIR, "loocv_all_predictions.csv"), row.names = FALSE)
saveRDS(conf_matrix, file.path(SUM_DIR, "loocv_conf_matrix.rds"))
saveRDS(overall_metrics, file.path(SUM_DIR, "loocv_overall_metrics.rds"))
saveRDS(models_list, file.path(SUM_DIR, "loocv_model_list.rds"))
saveRDS(roc_obj, file.path(SUM_DIR, "loocv_roc_curve.rds"))
saveRDS(imp_results, file.path(SUM_DIR, "loocv_imp_results.rds"))

print("LOOCV Complete")
print(overall_metrics)
print(conf_matrix)
