library(methylModes)

# Essential packages for Shiny app functionality
library(ggplot2)
library(plotly)
library(data.table)
library(dplyr)
library(DT)
library(parallelly)
library(foreach)
library(iterators)
library(shinyjs)
library(shinycssloaders)
library(shinybusy)
library(shinyWidgets)

# Load configuration values
source(system.file("shiny", "config.R", package = "methylModes"))
