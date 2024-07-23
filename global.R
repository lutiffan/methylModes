# General data manipulation
library(data.table)
library(dplyr)
library(tools) # file_ext()

# Shiny-specific
library(shiny)
library(shinyjs)
library(shinybusy)
library(shinycssloaders)
library(shinyWidgets)

# Tables and plots
library(ggplot2)
library(plotly)
library(DT)

# Parallelization
library(foreach)
library(iterators)
library(parallelly)
library(doParallel)

# Progress bar animation that handles parallel processes
library(doFuture)
library(progressr)

# Set up how future() calls are resolved
plan(multisession)
handlers(global = TRUE)
handlers("progress")

# Load default parameters
source(file = "standardParams.R", local = T)

# Load functions
source(file = "localMinMax.R", local = T)
source(file = "methylModes.R", local = T)
source(file = "fillPeakSummaryParallelProgressr.R", local = T)
# source("/home/lutiffan/peakDetectionScripts/fillPeakSummaryParallel.R")
source(file = "visualizeBetaPeaks.R", local = T)

# Global variables
KERNEL_TYPE = "gaussian"
BANDWIDTH_TYPE = NA
NUM_BREAKS = 500