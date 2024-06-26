# Shiny-specific
library(shiny)
library(shinyjs)
library(shinybusy)
library(shinycssloaders)

# Tables and plots
library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)

# Parallelization
library(foreach)
library(iterators)
library(parallelly)
library(doParallel)

# Load default parameters
source(file = "standardParams.R", local = T)

# Load functions
source(file = "localMinMax.R", local = T)
source(file = "methylModes.R", local = T)
source(file = "fillPeakSummaryParallel.R", local = T)
# source("/home/lutiffan/peakDetectionScripts/fillPeakSummaryParallel.R")
source(file = "visualizeBetaPeaks.R", local = T)

# Global variables
KERNEL_TYPE = "gaussian"
BANDWIDTH_TYPE = NA
NUM_BREAKS = 500