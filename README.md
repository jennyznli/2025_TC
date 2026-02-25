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
- `mir200c.R`: Infer miR-200c with CytoMethIC (1.6.0)
- `cell_deconvolution.R`: Cell-type deconvolution with EpiDISH (2.26.0)

### Clustering analysis
- `pca.R`: PCA on most variable CpGs
- `ped_clustering.R`, `adt_clustering.R`, `jnt_clustering.R`: Consensus clustering and generation of t-SNE embeddings
- `ped_tsne.R`, `adt_clustering.R`, `jnt_clustering.R`: Plot t-SNE visualizations (Fig. 1)

### Differential methylation
- `ped_diff_meth_inv.R`: Differential methylation on pediatric cohort by invasiveness
- `ped_diff_meth_cluster_leuko.R`, `ped_diff_meth_cluster_li.R`: Differential methylation on pediatric cohort by cluster
- `adt_diff_meth_inv.R`: Differential methylation on adult cohort by invasiveness
- `ped_diff_meth_analysis.R`, `adt_diff_meth_analysis.R`: Plot differential methylation analysis (Fig. 2)

### RNA analysis
- `rna.R`: RNA-seq analysis including differential expression, integrated RNA and methylation transcription factor binding site analysis, and GSEA. 
  
### Age analysis 
- `jnt_age.R`: Epigenetic clock inference with methylclock (1.16.0) and joint age analysis (Fig. 3). 

### Development of driver and invasiveness classifiers
- `loocv_invasiveness.R`, `loocv_driver.R`: Leave one out cross validation (LOOCV) on the reference cohort.
- `loocv_invasiveness_analysis.R`, `loocv_driver_analysis.R`: Evaluate and plot LOOCV performance (Fig. 4). 
- `fmodel_invasiveness.R`, `fmodel_driver.R`: Final model training on the validation cohort, testing on the validation cohort, evaluate performance, enrichment plots, and SHAP analysis (Fig. 4).  
- `fmodel_analysis.R`:  Plot cross model analysis figures (Fig. 4)
- `ped_stats.R`: Plot cohort statistical summary figures (Fig. 4). 

## Support

**Contact**: [jennyzli\@sas.upenn.edu](mailto:jennyzli@sas.upenn.edu)

