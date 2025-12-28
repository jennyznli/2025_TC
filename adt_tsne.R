# ============================================================
#   This script plots tSNE visualizations for the adult cohort.
# ============================================================
source("../config.R")

PLOT_SIZE <- list(width = 5.5, height = 4.5)

# ========================
# LOAD DATA
# ========================
ss <- readxl::read_excel(file.path("ss", ADT_META))

tsne_all <- read.csv(file.path("data", "adt_tsne_all_coords.csv"))

ss$Lymph_Node <- as.factor(ss$Lymph_Node)
ss$Chronological_Age <- as.numeric(ss$Chronological_Age)
ss$M <- as.factor(ss$M)
ss$N <- as.factor(ss$N)
ss$Probe_Success_Rate <- as.numeric(ss$Probe_Success_Rate)
ss$IC_EpiDISH <- as.numeric(ss$IC_EpiDISH)
ss$ERK_Score <- as.numeric(ss$ERK_Score)
ss$BRS_Score <- as.numeric(ss$BRS_Score)
ss$Differentiation_Score <- as.numeric(ss$Differentiation_Score)

overlap_cols <- intersect(names(tsne_primary), names(ss))
overlap_cols <- setdiff(overlap_cols, "IDAT")

tsne_all_clean <- tsne_all %>% select(-all_of(overlap_cols))
tsne_all <- inner_join(tsne_all_clean, ss, by = "IDAT")

# ========================
# PLOT CONFIGURATIONS
# ========================
plot_list <- list(
  list(name = "invasiveness", color_var = "Clinical_Invasiveness", color_values = invasiveness_colors),
  list(name = "driver_group", color_var = "Driver_Group", color_values = driver_colors),
  list(name = "sex", color_var = "Sex", color_values = sex_colors),
  list(name = "lymph_node", color_var = "Lymph_Node", color_values = lymph_node_colors),
  list(name = "probe_success_rate", color_var = "Probe_Success_Rate", continuous = TRUE),
  list(name = "leukocyte_fraction", color_var = "IC_EpiDISH", continuous = TRUE),
  list(name = "horvath_age", color_var = "Epigenetic_Age", continuous = TRUE),
  list(name = "N", color_var = "N", color_values = n_colors),
  list(name = "M", color_var = "M", color_values = m_colors),
  list(name = "ERK", color_var = "ERK_Score", continuous = TRUE),
  list(name = "BRS", color_var = "BRS_Score", continuous = TRUE),
  list(name = "Differentiation", color_var = "Differentiation_Score", continuous = TRUE)
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
  
  # alpha transparency for confidence
  if (!is.null(alpha_var)) {
    plot_data <- data[!is.na(data[[alpha_var]]), ]
    na_data <- data[is.na(data[[alpha_var]]), ]
    if (!is.null(shape_var)) {
      p <- p + geom_point(data = plot_data, aes(color = !!sym(color_var), shape = !!sym(shape_var), alpha = !!sym(alpha_var)), size = 2) +
        scale_shape_manual(values = shape_values) +
        scale_alpha_continuous(range = c(0.4, 1))
    } else {
      p <- p + geom_point(data = plot_data, aes(color = !!sym(color_var), alpha = !!sym(alpha_var)), size = 2, shape = 16) +
        scale_alpha_continuous(range = c(0.4, 1))
    }
    
    if (nrow(na_data) > 0) {
      if (!is.null(shape_var)) {
        p <- p + geom_point(data = na_data, aes(color = !!sym(color_var), shape = !!sym(shape_var)),
                            size = 2, alpha = 0.2)
      } else {
        p <- p + geom_point(data = na_data, aes(color = !!sym(color_var)),
                            size = 2, shape = 16, alpha = 0.2)
      }
    }
  } else {
    if (!is.null(shape_var)) {
      p <- p + geom_point(aes(color = !!sym(color_var), shape = !!sym(shape_var)), size = 2) +
        scale_shape_manual(values = shape_values)
    } else {
      p <- p + geom_point(aes(color = !!sym(color_var)), size = 2, shape = 16)
    }
  }
  
  # color scale
  if (continuous) {
    p <- p + scale_color_viridis_c()
  } else {
    p <- p + scale_color_manual(values = color_values)
  }
  
  # legend
  if (show_legend) {
    p <- p + theme(legend.position = "right", legend.text = element_text(size = 10), legend.key.size = unit(0.5, "cm"))
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

# ========================
# STANDARD PLOTS
# ========================
data_subsets <- list(
  list(name = "primary", data = tsne_all)
)

for (subset in data_subsets) {
  for (config in plot_list) {
    # without legend
    plot <- create_tsne_plot(
      data = subset$data,
      color_var = config$color_var,
      color_values = config$color_values,
      shape_var = config$shape_var,
      shape_values = config$shape_values,
      continuous = if (!is.null(config$continuous)) config$continuous else FALSE,
      show_legend = FALSE
    )
    
    filename <- file.path("figures", sprintf("adt_%s_tsne_%s.pdf", subset$name, config$name))
    ggsave(filename, plot, width = PLOT_SIZE$width, height = PLOT_SIZE$height, units = "in")
    
    # with legend
    plot_legend <- create_tsne_plot(
      data = subset$data,
      color_var = config$color_var,
      color_values = config$color_values,
      shape_var = config$shape_var,
      shape_values = config$shape_values,
      continuous = if (!is.null(config$continuous)) config$continuous else FALSE,
      show_legend = TRUE
    )
    
    filename_legend <- file.path("figures", sprintf("adt_%s_tsne_%s_legend.pdf", subset$name, config$name))
    ggsave(filename_legend, plot_legend, width = PLOT_SIZE$width, height = PLOT_SIZE$height, units = "in")
  }
}

