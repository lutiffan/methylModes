library(methylModes)

# # General data manipulation
# library(data.table)
# library(dplyr)
# library(tools) # file_ext()
# # library(qs) # qread()
# 
# # Shiny-specific
# library(shiny)
# library(shinyjs)
# library(shinybusy)
# library(shinycssloaders)
# library(shinyWidgets)
# library(shinythemes)
# 
# # Tables and plots
# library(ggplot2)
# library(plotly)
# library(DT)
# 
# # Parallelization
# library(foreach)
# library(iterators)
# library(parallelly)
# library(doParallel)
# 
# # Progress bar animation that handles parallel processes
# library(progress)
# library(doFuture)
# library(progressr)

# Set up how future() calls are resolved
# plan(multisession)
handlers(global = TRUE)
handlers("progress")

# Load default parameters from the package
# The parameters are now available through the package functions
# No need to source external files

# Global variables (these should be set through function parameters)
KERNEL_TYPE = "gaussian"
BANDWIDTH_TYPE = NA
NUM_BREAKS = 500
