# MASTER SCRIPT 

### PREPROCESSING ### 
# extract methylation betas from IDAT
extract_idats.R ped QCDPB FALSE
extract_idats.R adt QCDPB FALSE

# get QC stats
stats.R ped
stats.R adt

# remove sex probes, filter CG probes
qc_probes.R ped_betas_QCDPB.rds EPICv2
qc_probes.R adt_betas_QCDPB.rds HM450

# join adult and pediatric cohorts
jnt_probes.R

# mir200c inference 
mir200c.R ped_betas_QCDPB_prc.rds EPIC
mir200c.R adt_betas_QCDPB_pr.rds EPIC

# cell deconvolution 
cell_deconvolution.R  

# PCA with 30k most variable features
pca.R ped_betas_QCDPB_prc.rds 30000
pca.R ped_betas_QCDPB_pr.rds 30000
pca.R jnt_betas_QCDPB_pr.rds 30000

### FIG. 1 ### 
ped_stats.R 

# clustering analysis 
ped_clustering.R
adt_clustering.R
jnt_clustering.R

# tSNE visualization 
ped_tsne.R
adt_tsne.R
jnt_tsne.R

### FIG. 2 ### 
# differential methylation
ped_diff_meth_inv.R
ped_diff_meth_cluster_leuko.R
ped_diff_meth_cluster_li.R
adt_diff_meth_inv.R

# differential methylation analysis 
ped_diff_meth_analysis.R 
adt_diff_meth_analysis.R 

# RNA analysis 
rna.R 

### FIG. 3 ### 
# age analysis 
jnt_age.R 

### FIG. 4 ### 
# LOOCV models 
loocv_invasiveness.R
loocv_driver.R
loocv_invasiveness_analysis.R
loocv_driver_analysis.R 

# final models training, testing, and analysis
fmodel_invasiveness.R
fmodel_driver.R
fmodel_analysis.R 



