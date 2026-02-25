# DNA Methylation-based Risk Stratification and Classification of Pediatric Thyroid Carcinoma Invasiveness
This repository is for raw data processing and downstream data analysis & visualization in the "DNA Methylation-based Risk Stratification and Classification of Pediatric Thyroid Carcinoma Invasiveness" manuscript.


## Repository Structure

``` plaintext
project_root/
├── scripts/                            # Analysis scripts
├── R/                                  # Helper scripts
├── models/                             # Models 
├── config.R                             
└── README.md                              
```
## Setup 
Edit the config with the appropriate file names. All data analysis and visualization were performed using R (4.5.2) and R package SeSAMe (1.28.1).

## Model Usage
DNA methylation data was obtained using Infinium MethylationEPIC. Standard processing, including background correction, normalization, etc. must be performed to obtained methylation beta values (0-1). An example processing pipeline with SeSaME is shown below. If using EPICv2.0, you must collapse the probe names down to the prefixes before applying the classifier. There must be no NAs in the data. 
```
model <- readRDS("models/invasiveness_model.rds")
sel_probes <- rownames(model$importance) # CpG probes used during training
betas_sel  <- betas[sel_probes, ] # CpGs x samples
pred_class <- predict(model, t(betas_sel))
pred_prob  <- predict(model, t(betas_sel), type = "prob") # view probabilities
```
The invasiveness model will predict "High" or "Low". 

The driver model will predict one of the four classes: "Kinase Fusion," "BRAFV600E," "Ras-like," or "DICER1."

## Data Analysis
The scripts were written in a manner so that they are sequentially named depending on the processing. Usage details are found in the script headers. 

### Preprocessing
Extract beta values from IDAT files, perform QC, and optional masking.   
```
Rscript extract_idats.R ped QCDPB FALSE
Rscript extract_idats.R adt QCDPB FALSE
```

Remove CpGs on sex chromosomes, filter down to CG probes, and collapses EPICv2 data to prefixes. 
```
Rscript qc_probes.R ped_betas_QCDPB.rds EPICv2
Rscript qc_probes.R adt_betas_QCDPB.rds HM450
```

Join pediatric and adult cohorts by subsetting to common probes.
```
Rscript jnt_probes.R
```

### Inference
Infer miR-200c status with CytoMethIC
```
Rscript mir200c.R ped_betas_QCDPB_prc.rds EPIC
Rscript mir200c.R adt_betas_QCDPB_pr.rds EPIC
```

Cell-type deconvolution with EpiDISH 
```
Rscript cell_deconvolution.R
```

### Unsupervised clustering
Run PCA on top 30k most variable CpGs
```
Rscript pca.R ped_betas_QCDPB_prc.rds 30000
Rscript pca.R adt_betas_QCDPB_pr.rds 30000
Rscript pca.R jnt_betas_QCDPB_pr.rds 30000
```
Consensus clustering and generation of t-SNE embeddings
```
Rscript ped_clustering.R
Rscript adt_clustering.R
Rscript jnt_clustering.R

### Figure 1
```
Plotting of t-SNE visualizations
```
Rscript ped_tsne.R
Rscript adt_tsne.R
Rscript jnt_tsne.R
```

### Differential methylation
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

### RNA analysis
```
RNA-seq analysis (differential expression, GSEA)
```
Rscript rna.R
```

### Age analysis 
Epigenetic age analysis
```
Rscript jnt_age.R
```

### Random forest development 
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

## Support

**Contact**: [jennyzli\@sas.upenn.edu](mailto:jennyzli@sas.upenn.edu)

