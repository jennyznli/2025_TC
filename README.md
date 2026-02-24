# DNA Methylation-based Risk Stratification and Classification of Pediatric Thyroid Carcinoma Invasiveness

## Repository Structure

``` plaintext
project_root/
├── scripts/                            # Analysis scripts
├── R/                                  # Helper scripts
├── models/                             # Models 
├── README.md                           
└── config.R                                 
```

## Overview

### Preprocessing
Extract beta values from IDAT files
```
Rscript extract_idats.R ped QCDPB FALSE
Rscript extract_idats.R adt QCDPB FALSE
```

Compute QC statistics
```
Rscript qc_probes.R ped_betas_QCDPB.rds EPICv2
Rscript qc_probes.R adt_betas_QCDPB.rds HM450
```

Join pediatric and adult cohorts
```
Rscript jnt_probes.R
```

Infer miR-200c status with CytoMethIC
```
Rscript mir200c.R ped_betas_QCDPB_prc.rds EPIC
Rscript mir200c.R adt_betas_QCDPB_pr.rds EPIC
```

Cell-type deconvolution with EpiDISH 
```
Rscript cell_deconvolution.R
```

Run PCA on top 30k most variable CpGs
```
Rscript pca.R ped_betas_QCDPB_prc.rds 30000
Rscript pca.R adt_betas_QCDPB_pr.rds 30000
Rscript pca.R jnt_betas_QCDPB_pr.rds 30000
```

### Figure 1
Consensus clustering and generation of t-SNE embeddings
```
Rscript ped_clustering.R
Rscript adt_clustering.R
Rscript jnt_clustering.R
```
Plotting of t-SNE visualizations
```
Rscript ped_tsne.R
Rscript adt_tsne.R
Rscript jnt_tsne.R
```

### Figure 2
Differential methylation on pediatric cohort (invasiveness, cluster)
```
Rscript ped_diff_meth_inv.R
Rscript ped_diff_meth_cluster_leuko.R
Rscript ped_diff_meth_cluster_li.R
```

Differential methylation on adult cohort (invasiveness)
```
Rscript adt_diff_meth_inv.R
```

Plotting and analysis of differential methylation
```
Rscript ped_diff_meth_analysis.R
Rscript adt_diff_meth_analysis.R

```
RNA-seq analysis (differential expression, GSEA)
```
Rscript rna.R
```

### Figure 3
Epigenetic age analysis
```
Rscript jnt_age.R
```

### Figure 4
Leave one out cross validation (LOOCV) on the reference cohort
```
Rscript loocv_invasiveness.R
Rscript loocv_driver.R
```

Plotting LOOCV figures
```
Rscript loocv_invasiveness_analysis.R
Rscript loocv_driver_analysis.R
```

Final model development with the reference cohort and testing on the validation cohort
```
Rscript fmodel_invasiveness.R
Rscript fmodel_driver.R
Rscript fmodel_analysis.R
```

Generation of cohort figures
```
Rscript ped_stats.R
```

### Model Usage
Download and load the model. QC your methylation betas, and ensure that you have no NAs and subset to the feature probes used by the model. Apply the model like so
```
model     <- readRDS("models/invasiveness_model.rds")
sel_probes <- model$finalModel$xNames   # probes used during training

betas_sel  <- betas[sel_probes, ]       # CpGs x samples

pred_class <- predict(model, t(betas_sel))
pred_prob  <- predict(model, t(betas_sel), type = "prob")

predictions <- data.frame(
  IDAT            = colnames(betas_sel),
  Predicted_Class = pred_class,
  pred_prob,
  Max_Probability = apply(pred_prob, 1, max)
)

write.csv(predictions, "invasiveness_predictions.csv", row.names = FALSE)
```

## Support

**Contact**: [jennyzli\@sas.upenn.edu](mailto:jennyzli@sas.upenn.edu)

