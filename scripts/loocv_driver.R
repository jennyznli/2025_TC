# ============================================================
#   This script performs leave-one-out cross validation (LOOCV)
#   to develop and validate a random forest classifier for
#   predicting driver group in thyroid cancer samples
#   using DNA methylation data.
# ============================================================
N_WORKERS <- 20
set.seed(123)

source("../config.R")

CL_DIR <- file.path(BASE_DIR, "data", "loocv_driver_model")
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
    n_batches <- ceiling(n_probes/batch_size)

    fold_models <- list()
    fold_importance <- data.frame()

    for(s in 0:(n_batches-1)) {
        # Calculate indices for current batch
        start_idx <- s * batch_size + 1
        end_idx <- min((s + 1) * batch_size, n_probes)
        sbetas <- betas[start_idx:end_idx, , drop=FALSE]
        print(paste("Batch dimensions:", dim(sbetas)[1], "x", dim(sbetas)[2]))

        model <- randomForest(
            x = t(sbetas[, train_idx]),
            y = as.factor(labels[train_idx]),
            ntree = n_trees,
            importance = TRUE
        )

        fold_models[[s+1]] <- model

        importance <- data.frame(
            Fold = if(is.null(fold)) NA else fold,
            Model = s,
            Feature = rownames(model$importance),
            MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"]
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

    # Model training
    model <- randomForest(
        x = t(betas[sel_probes, train_idx]),
        y = labels[train_idx],
        ntree = 500,
        importance = TRUE
    )

    # importance
    importance <- data.frame(
        Feature = rownames(model$importance),
        MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"],
        MeanDecreaseGini = model$importance[,"MeanDecreaseGini"]
    )

    # Add class-specific importance
    for (level in levels(labels)) {
        if(level %in% colnames(model$importance)) {
            importance[,level] <- model$importance[,level]
        }
    }

    importance$Fold <- i

    # Predictions
    pred_probs <- predict(model, t(betas[sel_probes, valid_idx, drop=FALSE]), type = "prob")
    pred_class <- predict(model, t(betas[sel_probes, valid_idx, drop=FALSE]))

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
        preds[,col_name] <- pred_probs[,j]
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
# EXECUTION
# ========================
ss <- read_excel(file.path("ss", PED_META)) %>%
    dplyr::filter(Lymph_Node == "F", Driver_Group != "Indeterminate")
betas <- readRDS(file.path("data", PED_BETAS))

valid_ids <- ss$IDAT[ss$IDAT %in% colnames(betas)]
ss <- ss[match(valid_ids, ss$IDAT), ]
betas <- betas[, valid_ids]

labels <- as.factor(ss$Driver_Group)

# Create sample indices for LOOCV
all_samples <- 1:length(labels)
print(paste("Total samples for LOOCV:", length(all_samples)))

print("Starting LOOCV")
results <- bplapply(1:length(all_samples), function(i) process_sample(i, all_samples, betas, labels), BPPARAM = MulticoreParam(N_WORKERS))
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

overall_metrics <- data.frame(
    Accuracy = accuracy,
    Overall_Kappa = conf_matrix$overall["Kappa"]
)

# multiclass ROC curves (one-vs-all)
roc_list <- list()
auc_values <- numeric(length(levels(labels)))

for (i in 1:length(levels(labels))) {
    level <- levels(labels)[i]
    safe_level_name <- gsub(" ", ".", gsub("-", ".", level))
    prob_col <- paste0("Prob_", safe_level_name)

    # verify the column exists
    if (!prob_col %in% colnames(pred_results)) {
        warning(paste("Column", prob_col, "not found in pred_results for level", level))
        print(paste("Available columns:", paste(colnames(pred_results), collapse=", ")))
        auc_values[i] <- NA
        next
    }

    # binary outcome (this class vs all others)
    binary_outcome <- ifelse(pred_results$True_Label == level, 1, 0)

    # Diagnostic information
    cat("Processing class:", level, "\n")
    cat("Using probability column:", prob_col, "\n")
    cat("Number of positive samples:", sum(binary_outcome), "\n")
    cat("Number of negative samples:", length(binary_outcome) - sum(binary_outcome), "\n")

    # Calculate ROC if we have at least one sample of this class
    if (sum(binary_outcome) > 0 && sum(binary_outcome) < length(binary_outcome)) {
        valid_indices <- !is.na(binary_outcome) & !is.na(pred_results[[prob_col]])
        if (sum(valid_indices) < 2) {
            warning("Not enough valid data points for ROC curve")
            auc_values[i] <- NA
            next
        }

        tryCatch({
            roc_obj <- roc(binary_outcome, pred_results[[prob_col]],
                           quiet = TRUE, direction = "<")  # Ensures proper direction
            roc_list[[level]] <- roc_obj
            auc_values[i] <- as.numeric(pROC::auc(roc_obj))
            coords <- coords(roc_obj, "all")
        }, error = function(e) {
            warning(paste("Error calculating ROC for class", level, ":", e$message))
            auc_values[i] <- NA
        })
    } else {
        warning(paste("Class", level, "does not have both positive and negative samples"))
        auc_values[i] <- NA
    }
}

# Add class AUC to overall metrics
for (i in 1:length(levels(labels))) {
    level <- levels(labels)[i]
    overall_metrics[[paste0("AUC_", level)]] <- auc_values[i]
}

# Calculate mean AUC across classes
overall_metrics$Mean_AUC <- mean(auc_values, na.rm = TRUE)

class_metrics <- as.data.frame(conf_matrix$byClass)

# Process per-class metrics using one-vs-all approach
per_class_metrics <- data.frame()

for (level in levels(labels)) {
    # Calculate metrics for this class vs all others
    binary_true <- ifelse(pred_results$True_Label == level, "Positive", "Negative")
    binary_pred <- ifelse(pred_results$Predicted == level, "Positive", "Negative")

    cat("Processing one-vs-all for class:", level, "\n")
    cat("Positive true cases:", sum(binary_true == "Positive"), "\n")
    cat("Positive predicted cases:", sum(binary_pred == "Positive"), "\n")

    tryCatch({
        binary_conf <- confusionMatrix(
            data = factor(binary_pred, levels = c("Positive", "Negative")),
            reference = factor(binary_true, levels = c("Positive", "Negative"))
        )
        class_row <- data.frame(
            Class = level,
            Sensitivity = binary_conf$byClass["Sensitivity"],
            Specificity = binary_conf$byClass["Specificity"],
            Precision = binary_conf$byClass["Pos Pred Value"],
            Recall = binary_conf$byClass["Sensitivity"],
            F1 = binary_conf$byClass["F1"],
            BalancedAccuracy = binary_conf$byClass["Balanced Accuracy"]
        )

        per_class_metrics <- rbind(per_class_metrics, class_row)
    }, error = function(e) {
        warning(paste("Error in confusion matrix for class", level, ":", e$message))
        class_row <- data.frame(
            Class = level,
            Sensitivity = NA,
            Specificity = NA,
            Precision = NA,
            Recall = NA,
            F1 = NA,
            BalancedAccuracy = NA
        )
        per_class_metrics <- rbind(per_class_metrics, class_row)
    })
}

macro_metrics <- data.frame(
    Class = "Macro-Average",
    Sensitivity = mean(per_class_metrics$Sensitivity, na.rm = TRUE),
    Specificity = mean(per_class_metrics$Specificity, na.rm = TRUE),
    Precision = mean(per_class_metrics$Precision, na.rm = TRUE),
    Recall = mean(per_class_metrics$Recall, na.rm = TRUE),
    F1 = mean(per_class_metrics$F1, na.rm = TRUE),
    BalancedAccuracy = mean(per_class_metrics$BalancedAccuracy, na.rm = TRUE)
)
per_class_metrics <- rbind(per_class_metrics, macro_metrics)

overall_metrics$Macro_Sensitivity <- macro_metrics$Sensitivity
overall_metrics$Macro_Specificity <- macro_metrics$Specificity
overall_metrics$Macro_Precision <- macro_metrics$Precision
overall_metrics$Macro_Recall <- macro_metrics$Recall
overall_metrics$Macro_F1 <- macro_metrics$F1
overall_metrics$Macro_BalancedAccuracy <- macro_metrics$BalancedAccuracy

# ========================
# SAVE FINAL RESULTS
# ========================
saveRDS(per_class_metrics, file.path(SUM_DIR, "loocv_per_class_metrics.rds"))
write.csv(pred_results, file.path(SUM_DIR, "loocv_all_predictions.csv"), row.names = FALSE)
saveRDS(conf_matrix, file.path(SUM_DIR, "loocv_conf_matrix.rds"))
saveRDS(overall_metrics, file.path(SUM_DIR, "loocv_overall_metrics.rds"))
saveRDS(models_list, file.path(SUM_DIR, "loocv_model_list.rds"))
saveRDS(roc_list, file.path(SUM_DIR, "loocv_roc_curves.rds"))
saveRDS(imp_results, file.path(SUM_DIR, "loocv_imp_results.rds"))

print("LOOCV Complete")
print(overall_metrics)
print(conf_matrix)
write.csv(overall_metrics, file.path(SUM_DIR, "loocv_overall_metrics.csv"), row.names = FALSE)
