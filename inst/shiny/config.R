#' Shiny App Configuration
#' 
#' This file contains configuration values for the methylModes Shiny application.
#' These values can be modified to customize the app behavior.

# File upload settings
MAX_FILE_SIZE <- 100 * 1024^3  # 100 GB

# Plot settings
DEFAULT_HISTOGRAM_BINS <- 50
DEFAULT_DENSITY_ADJUST <- 1.5
DEFAULT_PUSH_TO_ZERO <- 1e-6

# Analysis thresholds
DEFAULT_VARIANCE_THRESHOLD <- 0.01
DEFAULT_HYPO_THRESHOLD <- 0.3
DEFAULT_HYPER_THRESHOLD <- 0.7

# Fixed algorithm parameters (these are not user-configurable in the UI)
FIXED_KERNEL_TYPE <- "gaussian"
FIXED_BANDWIDTH_TYPE <- NA
FIXED_NUM_BREAKS <- 500

# Color schemes for plots
ISLAND_COLORS <- c(
  "Island" = "#FDE333",
  "N_Shore" = "#C1DE35", 
  "S_Shore" = "#00C376",
  "N_Shelf" = "#008498",
  "S_Shelf" = "#006892",
  "OpenSea" = "#363D7C"
)

DENSITY_PLOT_COLORS <- c(
  "Density Estimate" = "#00588B",
  "Maxima" = "#B2DC3C",
  "Minima (peak boundaries)" = "#009B95"
) 