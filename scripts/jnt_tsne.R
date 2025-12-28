# ============================================================
#   This script plots tSNE visualizations for the joint cohort.
# ============================================================
source("../config.R")

PLOT_SIZE <- list(width = 5.5, height = 4.5)

# ========================
# LOAD DATA
# ========================
ss <- readxl::read_excel(file.path("ss", JNT_META))
tsne_primary_ref <- read.csv(file.path("data", "jnt_tsne_primary_ref_coords.csv"))

ss$Lymph_Node <- as.factor(ss$Lymph_Node)
ss$Chronological_Age <- as.numeric(ss$Chronological_Age)
ss$M <- as.factor(ss$M)
ss$N <- as.factor(ss$N)
ss$Probe_Success_Rate <- as.numeric(ss$Probe_Success_Rate)
ss$IC_EpiDISH <- as.numeric(ss$IC_EpiDISH)
ss$Source <- as.factor(ss$Source)
ss$Reference <- as.factor(ss$Reference)

tsne_primary_ref <- dplyr::left_join(tsne_primary_ref, ss, by = "IDAT")

# ========================
# PLOT CONFIGURATIONS
# ========================
plot_list <- list(
  list(name = "invasiveness", color_var = "Clinical_Invasiveness", color_values = invasiveness_colors),
  list(name = "driver_group", color_var = "Driver_Group", color_values = driver_colors),
  list(name = "sex", color_var = "Sex", color_values = sex_colors),
  list(name = "cluster", color_var = "Consensus_Cluster", color_values = cluster_colors),
  list(name = "lymph_node", color_var = "Lymph_Node", color_values = lymph_node_colors),
  list(name = "probe_success_rate", color_var = "Probe_Success_Rate", continuous = TRUE),
  list(name = "leukocyte_fraction", color_var = "IC_EpiDISH", continuous = TRUE),
  list(name = "horvath_age", color_var = "Epigenetic_Age", continuous = TRUE),
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
create_tsne_plot <- function(data, color_var, color_values = NULL, continuous = FALSE, 
                             show_legend = FALSE) {
  
  p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2)) +
    coord_fixed() +
    common_theme
  
  p <- p + geom_point(aes(color = !!sym(color_var), shape = Source), size = 2) +
    scale_shape_manual(values = source_shapes)
  
  if (continuous) {
    p <- p + scale_color_viridis_c()
  } else {
    p <- p + scale_color_manual(values = color_values)
  }
  
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
  list(name = "primary_ref", data = tsne_primary_ref)
)

for (subset in data_subsets) {
  for (config in plot_list) {
    # without legend
    plot <- create_tsne_plot(
      data = subset$data,
      color_var = config$color_var,
      color_values = config$color_values,
      continuous = if (!is.null(config$continuous)) config$continuous else FALSE,
      show_legend = FALSE
    )
    filename <- file.path("figures", sprintf("jnt_%s_tsne_%s.pdf", subset$name, config$name))
    ggsave(filename, plot, width = PLOT_SIZE$width, height = PLOT_SIZE$height, units = "in")
    
    # with legend
    plot_legend <- create_tsne_plot(
      data = subset$data,
      color_var = config$color_var,
      color_values = config$color_values,
      continuous = if (!is.null(config$continuous)) config$continuous else FALSE,
      show_legend = TRUE
    )
    
    filename_legend <- file.path("figures", sprintf("jnt_%s_tsne_%s_legend.pdf", subset$name, config$name))
    ggsave(filename_legend, plot_legend, width = PLOT_SIZE$width, height = PLOT_SIZE$height, units = "in")
  }
}
