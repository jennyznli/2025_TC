# ============================================================
# t-SNE PLOTS
# ============================================================
source("../config.R")

PLOT_SIZE <- list(width = 5.5, height = 4.5)

# ========================
# MERGE BEFOREHAND...
# ========================
# # RNA expression for ZEB1/2
# ss <- read_excel(file.path("ss", PED_META))
# rna <- read.csv(file.path("data", "20251209_rna_log2cpmfiltered.csv"), row.names = 1)
# rna <- as.data.frame(t(rna)) %>% select("ZEB1", "ZEB2")
# rna$Sample_ID <- rownames(rna)
# 
# # MIR200C
# mir200c <- as.data.frame(readRDS(file.path("data", "ped_betas_QCDPB_prc_mir200c.rds")))
# mir200c$IDAT <- rownames(mir200c)
# 
# # join everything... 
# ss <- left_join(ss, mir200c, by="IDAT") 
# ss <- left_join(ss, rna, by="Sample_ID")
# write.csv(ss, file.path("ss", "test.csv"))

# ========================
# LOAD DATA
# ========================
ss <- read_excel(file.path("ss", PED_META)) %>% filter(LN_Include_In_Analysis == 1)

tsne_primary_1 <- read.csv(file.path("data", "ped_tsne_primary_1_coords.csv"))
tsne_primary_12 <- read.csv(file.path("data", "ped_tsne_primary_12_coords.csv"))
tsne_all_1 <- read.csv(file.path("data", "ped_tsne_all_1_coords.csv"))
tsne_all_12 <- read.csv(file.path("data", "ped_tsne_all_12_coords.csv"))
ss$Lymph_Node <- as.factor(ss$Lymph_Node)
ss$Chronological_Age <- as.numeric(ss$Chronological_Age)
ss$M <- as.factor(ss$M)
ss$N <- as.factor(ss$N)
ss$'T' <- as.factor(ss$'T')
ss$Tissue_Type <- as.factor(ss$Tissue_Type)
ss$Batch <- as.factor(ss$Batch)
ss$Probe_Success_Rate <- as.numeric(ss$Probe_Success_Rate)
ss$IC_EpiDISH <- as.numeric(ss$IC_EpiDISH)
ss$ZEB1 <- as.numeric(ss$ZEB1)
ss$ZEB2 <- as.numeric(ss$ZEB2)
ss$MIR200C <- as.numeric(ss$MIR200C)

ss$Fmodel_Invasiveness_Predicted <- as.factor(ss$Fmodel_Invasiveness_Predicted)
ss$Fmodel_Invasiveness_Confidence <- as.numeric(ss$Fmodel_Invasiveness_Confidence)
ss$Fmodel_Invasiveness_Accuracy <- as.factor(ss$Fmodel_Invasiveness_Accuracy)
ss$Fmodel_Driver_Predicted <- as.factor(ss$Fmodel_Driver_Predicted)
ss$Fmodel_Driver_Confidence <- as.numeric(ss$Fmodel_Driver_Confidence)
ss$Fmodel_Driver_Accuracy <- as.factor(ss$Fmodel_Driver_Accuracy)

ss$LOOCV_Invasiveness_Prediction <- as.factor(ss$LOOCV_Invasiveness_Prediction)
ss$LOOCV_Invasiveness_Confidence <- as.numeric(ss$LOOCV_Invasiveness_Confidence)
ss$LOOCV_Invasiveness_Accuracy <- as.factor(ss$LOOCV_Invasiveness_Accuracy)
ss$LOOCV_Driver_Prediction <- as.factor(ss$LOOCV_Driver_Prediction)
ss$LOOCV_Driver_Confidence <- as.numeric(ss$LOOCV_Driver_Confidence)
ss$LOOCV_Driver_Accuracy <- as.factor(ss$LOOCV_Driver_Accuracy)

tsne_primary_1 <- dplyr::left_join(tsne_primary_1, ss, by = "IDAT")
tsne_primary_12 <- dplyr::left_join(tsne_primary_12, ss, by = "IDAT")
tsne_all_1 <- dplyr::left_join(tsne_all_1, ss, by = "IDAT")
tsne_all_12 <- dplyr::left_join(tsne_all_12, ss, by = "IDAT")
            
# ========================
# PLOT CONFIGURATIONS
# ========================
plot_list <- list(
  list(name = "invasiveness", color_var = "Clinical_Invasiveness", color_values = invasiveness_colors),
  list(name = "driver_group", color_var = "Driver_Group", color_values = driver_colors),
  list(name = "sex", color_var = "Sex", color_values = sex_colors),
  list(name = "cluster", color_var = "CC_Cluster", color_values = cluster_colors),
  list(name = "batch", color_var = "Batch", color_values = batch_colors),
  list(name = "tissue_type", color_var = "Tissue_Type", color_values = tissue_colors),
  list(name = "lymph_node", color_var = "Lymph_Node", color_values = lymph_node_colors),
  list(name = "probe_success_rate", color_var = "Probe_Success_Rate", continuous = TRUE),
  list(name = "leukocyte_fraction", color_var = "IC_EpiDISH", continuous = TRUE),
  list(name = "horvath_age", color_var = "Horvath_MethylClock", continuous = TRUE),
  list(name = "age", color_var = "Chronological_Age", continuous = TRUE),
  list(name = "zeb1", color_var = "ZEB1", continuous = TRUE),
  list(name = "zeb2", color_var = "ZEB2", continuous = TRUE),
  list(name = "mir200c", color_var = "MIR200C", continuous = TRUE),
  list(name = "T", color_var = "T", color_values = t_colors),
  list(name = "N", color_var = "N", color_values = n_colors),
  list(name = "M", color_var = "M", color_values = m_colors)
)

common_theme <- theme(
  legend.title = element_blank(),
  panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray"),
  aspect.ratio = 1,
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
)

# ========================
# PLOTTING FUNCTIONS
# ========================
create_tsne_plot <- function(data, color_var, color_values = NULL, shape_var = NULL,
                             shape_values = NULL, continuous = FALSE, alpha_var = NULL,
                             show_legend = FALSE) {
  
  p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2)) +
    coord_fixed() +
    common_theme
  
  # Handle alpha variable for confidence visualization
  if (!is.null(alpha_var)) {
    plot_data <- data[!is.na(data[[alpha_var]]), ]
    na_data <- data[is.na(data[[alpha_var]]), ]
    
    if (!is.null(shape_var)) {
      p <- p + geom_point(data = plot_data,
                          aes(color = !!sym(color_var),
                              shape = !!sym(shape_var),
                              alpha = !!sym(alpha_var)),
                          size = 2) +
        scale_shape_manual(values = shape_values) +
        scale_alpha_continuous(range = c(0.4, 1))
    } else {
      p <- p + geom_point(data = plot_data,
                          aes(color = !!sym(color_var),
                              alpha = !!sym(alpha_var)),
                          size = 2, shape = 16) +
        scale_alpha_continuous(range = c(0.4, 1))
    }
    
    if (nrow(na_data) > 0) {
      if (!is.null(shape_var)) {
        p <- p + geom_point(data = na_data,
                            aes(color = !!sym(color_var),
                                shape = !!sym(shape_var)),
                            size = 2, alpha = 0.2)
      } else {
        p <- p + geom_point(data = na_data,
                            aes(color = !!sym(color_var)),
                            size = 2, shape = 16, alpha = 0.2)
      }
    }
  } else {
    if (!is.null(shape_var)) {
      p <- p + geom_point(aes(color = !!sym(color_var),
                              shape = !!sym(shape_var)),
                          size = 2) +
        scale_shape_manual(values = shape_values)
    } else {
      p <- p + geom_point(aes(color = !!sym(color_var)),
                          size = 2, shape = 16)
    }
  }
  
  # Apply color scale
  if (continuous) {
    p <- p + scale_color_viridis_c()
  } else {
    p <- p + scale_color_manual(values = color_values)
  }
  
  # Legend position
  if (show_legend) {
    p <- p + theme(legend.position = "right",
                   legend.text = element_text(size = 10),
                   legend.key.size = unit(0.5, "cm"))
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

# ========================
# STANDARD PLOTS
# ========================
data_subsets <- list(
  list(name = "primary_1", data = tsne_primary_1)
  # list(name = "primary_12", data = tsne_primary_12),
  # list(name = "all_1", data = tsne_all_1),
  # list(name = "all_12", data = tsne_all_12)
)

# Loop through each data subset
for (subset in data_subsets) {
  cat("Creating plots for:", subset$name, "\n")
  
  # Loop through each plot configuration
  for (config in plot_list) {
    # Plot without legend
    plot <- create_tsne_plot(
      data = subset$data,
      color_var = config$color_var,
      color_values = config$color_values,
      shape_var = config$shape_var,
      shape_values = config$shape_values,
      continuous = if (!is.null(config$continuous)) config$continuous else FALSE,
      show_legend = FALSE
    )
    
    filename <- file.path("figures", sprintf("ped_%s_tsne_%s.pdf", subset$name, config$name))
    ggsave(filename, plot, width = PLOT_SIZE$width, height = PLOT_SIZE$height, units = "in")
    
    # Plot with legend
    plot_legend <- create_tsne_plot(
      data = subset$data,
      color_var = config$color_var,
      color_values = config$color_values,
      shape_var = config$shape_var,
      shape_values = config$shape_values,
      continuous = if (!is.null(config$continuous)) config$continuous else FALSE,
      show_legend = TRUE
    )
    
    filename_legend <- file.path("figures", sprintf("ped_%s_tsne_%s_legend.pdf", subset$name, config$name))
    ggsave(filename_legend, plot_legend, width = PLOT_SIZE$width, height = PLOT_SIZE$height, units = "in")
  }
}

# ========================
# CONFIDENCE PLOTS - FMODEL 
# ========================

data_subsets <- list(
  list(name = "all_12", data = tsne_all_12)
)

confidence_plots <- list(
  list(
    name = "invasiveness_confidence",
    color_var = "Clinical_Invasiveness",
    color_values = invasiveness_colors,
    shape_var = "Fmodel_Invasiveness_Accuracy",
    shape_values = accuracy_shapes,
    alpha_var = "Fmodel_Invasiveness_Confidence"
  ),
  list(
    name = "drivergroup_confidence",
    color_var = "Driver_Group",
    color_values = driver_colors,
    shape_var = "Fmodel_Driver_Accuracy",
    shape_values = accuracy_shapes,
    alpha_var = "Fmodel_Driver_Confidence"
  )
)

create_confidence_plot <- function(data, color_var, color_values, shape_var, 
                                   shape_values, alpha_var, show_legend = FALSE) {
  
  # Separate data with and without confidence scores
  has_confidence <- data[!is.na(data[[alpha_var]]), ]
  no_confidence <- data[is.na(data[[alpha_var]]), ]
  
  p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2)) +
    coord_fixed() +
    common_theme
  
  # Plot samples WITHOUT confidence scores first (hollow circles, no shape variation)
  if (nrow(no_confidence) > 0) {
    p <- p + geom_point(data = no_confidence,
                        aes(color = !!sym(color_var), alpha = 0.7),
                        size = 2,
                        shape = 1)  # 1 = hollow circle
  }
  
  # Plot samples WITH confidence scores (filled shapes with alpha based on confidence)
  if (nrow(has_confidence) > 0) {
    p <- p + geom_point(data = has_confidence,
                        aes(color = !!sym(color_var),
                            shape = !!sym(shape_var),
                            alpha = !!sym(alpha_var)),
                        size = 2) +
      scale_shape_manual(values = shape_values) +
      scale_alpha_continuous(range = c(0.4, 1))
  }
  
  # Apply color scale
  p <- p + scale_color_manual(values = color_values)
  
  # Legend position
  if (show_legend) {
    p <- p + theme(legend.position = "right",
                   legend.text = element_text(size = 10),
                   legend.key.size = unit(0.5, "cm"))
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

# Loop through each data subset for confidence plots
for (subset in data_subsets) {
  cat("Creating confidence plots for:", subset$name, "\n")
  
  for (config in confidence_plots) {
    plot <- create_confidence_plot(
      data = subset$data,
      color_var = config$color_var,
      color_values = config$color_values,
      shape_var = config$shape_var,
      shape_values = config$shape_values,
      alpha_var = config$alpha_var,
      show_legend = FALSE
    )
    
    filename <- file.path("figures", sprintf("ped_%s_tsne_%s.pdf", subset$name, config$name))
    ggsave(filename, plot, width = PLOT_SIZE$width, height = PLOT_SIZE$height, units = "in")
  }
}
# ========================
# LOOCV CONFIDENCE PLOTS
# ========================
loocv_plot_configs <- list(
  list(
    name = "invasiveness_confidence",
    color_var = "Clinical_Invasiveness",
    color_values = invasiveness_colors,
    shape_var = "LOOCV_Invasiveness_Accuracy",
    shape_values = accuracy_shapes,
    continuous = FALSE,
    alpha_var = "LOOCV_Invasiveness_Confidence"
  ),
  list(
    name = "drivergroup_confidence",
    color_var = "Driver_Group",
    color_values = driver_colors,
    shape_var = "LOOCV_Driver_Accuracy",
    shape_values = accuracy_shapes,
    continuous = FALSE,
    alpha_var = "LOOCV_Driver_Confidence"
  )
)

for (config in loocv_plot_configs) {
  tryCatch({
    # Create plot without legend
    plot <- create_confidence_plot(
      data = tsne_primary_1,
      color_var = config$color_var,
      color_values = config$color_values,
      shape_var = config$shape_var,
      shape_values = config$shape_values,
      alpha_var = config$alpha_var,
      show_legend = FALSE
    )
    
    filename <- file.path("figures", sprintf("ped_primary_1_tsne_loocv_%s.pdf", config$name))
    ggsave(filename, plot, width = PLOT_SIZE$width, height = PLOT_SIZE$height, units = "in")
    
    # Create plot with legend
    plot_legend <- create_confidence_plot(
      data = tsne_primary_1,
      color_var = config$color_var,
      color_values = config$color_values,
      shape_var = config$shape_var,
      shape_values = config$shape_values,
      alpha_var = config$alpha_var,
      show_legend = TRUE
    )
    
    filename_legend <- file.path("figures", sprintf("ped_primary_1_tsne_loocv_%s_legend.pdf", config$name))
    ggsave(filename_legend, plot_legend, width = PLOT_SIZE$width, height = PLOT_SIZE$height, units = "in")
    
    cat("Successfully created LOOCV plot:", config$name, "\n")
  }, error = function(e) {
    cat("Error creating LOOCV plot for", config$name, ":", conditionMessage(e), "\n")
  })
}

# ========================
# PAIRED LYMPH NODE PLOT
# ========================
lymph_nodes <- tsne_all_12 %>%
  filter(Lymph_Node == "T") 

pairs_data <- data.frame()
pair_id <- 1

for (i in 1:nrow(lymph_nodes)) {
  ln <- lymph_nodes[i, ]
  if (is.null(ln$Paired_Primary) || is.na(ln$Paired_Primary) || ln$Paired_Primary == "") {
    next
  }
  
  primary <- tsne_all_12 %>%
    filter(Sample_ID == ln$Paired_Primary)
  
  if (nrow(primary) > 0) {
    ln_row <- ln %>%
      mutate(pair_id = paste0("Pair ", pair_id),
             sample_type = "Lymph Node")
    
    primary_row <- primary %>%
      mutate(pair_id = paste0("Pair ", pair_id),
             sample_type = "Primary Tumor")
    
    pairs_data <- rbind(pairs_data, ln_row, primary_row)
    pair_id <- pair_id + 1
  }
}

# Add pair ID for coloring - use the primary tumor ID for consistency
pairs_data <- pairs_data %>%
  mutate(Pair_ID = ifelse(sample_type == "Lymph Node", 
                          Paired_Primary, 
                          Sample_ID))

# Create connection lines - one line per LN-primary pair
connections <- data.frame()
for (current_pair_id in unique(pairs_data$pair_id)) {
  pair_samples <- pairs_data %>% filter(pair_id == current_pair_id)
  
  # Get lymph node and primary for this specific pair_id
  ln <- pair_samples %>% filter(sample_type == "Lymph Node")
  primary <- pair_samples %>% filter(sample_type == "Primary Tumor")
  
  # Make sure we have both a lymph node and primary
  if (nrow(ln) > 0 && nrow(primary) > 0) {
    connections <- rbind(connections, data.frame(
      pair_id = current_pair_id,
      Pair_ID = ln$Pair_ID,
      x1 = ln$tSNE1,
      y1 = ln$tSNE2,
      x2 = primary$tSNE1,
      y2 = primary$tSNE2
    ))
  }
}

# Remove duplicate primary tumors from pairs_data for plotting
pairs_data <- pairs_data %>%
  distinct(Sample_ID, .keep_all = TRUE)

# Create color mapping for pairs
n_pairs <- 15
pair_color_mapping <- setNames(
  pair_colors[1:n_pairs],
  unique(pairs_data$Pair_ID)
)

# Mark which samples are paired
tsne_all_12 <- tsne_all_12 %>%
  mutate(is_paired = Sample_ID %in% pairs_data$Sample_ID)

# Create plot
p_paired <- ggplot() +
  # Background points for non-paired samples
  geom_point(data = tsne_all_12 %>% filter(!is_paired),
             aes(x = tSNE1, y = tSNE2),
             color = "grey80", alpha = 0.4, size = 1.5) +
  # Connection lines
  geom_segment(data = connections,
               aes(x = x1, y = y1, xend = x2, yend = y2, color = Pair_ID),
               linetype = "dashed", linewidth = 0.7) +
  # Points for paired samples
  geom_point(data = pairs_data,
             aes(x = tSNE1, y = tSNE2, color = Pair_ID,
                 shape = sample_type, size = sample_type),
             alpha = 0.8) +
  scale_color_manual(values = pair_color_mapping, name = "Sample ID") +
  scale_shape_manual(values = c("Lymph Node" = 17, "Primary Tumor" = 1),
                     name = "Sample Type") +
  scale_size_manual(values = c("Lymph Node" = 3, "Primary Tumor" = 3.5),
                    name = "Sample Type") +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3))
  ) +
  coord_fixed() +
  common_theme +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

output_file <- file.path("figures", "ped_tsne_paired_ln.pdf")
ggsave(output_file, p_paired, width = 7, height = 5, units = "in")






