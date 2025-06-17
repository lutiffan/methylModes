#' Launch the methylModes Shiny Application
#' 
#' This function launches the interactive Shiny web application for analyzing
#' DNA methylation data using the methylModes algorithm. The app provides a
#' user-friendly interface for uploading data, running analyses, and visualizing
#' results.
#' 
#' @param launch.browser Logical. If TRUE (default), the app will open in the 
#'   default web browser.
#' @param port Integer. The port on which the app should run (default: 3838).
#' @param host Character. The host address (default: "127.0.0.1").
#' 
#' @return Launches a Shiny application in the user's default web browser.
#' 
#' @details
#' The Shiny app provides the following features:
#' \itemize{
#'   \item Upload and validate methylation data files
#'   \item Interactive data exploration and visualization
#'   \item Run methylModes analysis on individual probes or batch processing
#'   \item View and download results
#'   \item Customize analysis parameters
#' }
#' 
#' @examples
#' \dontrun{
#' # Launch the app
#' runMethylModesApp()
#' 
#' # Launch on a specific port
#' runMethylModesApp(port = 8080)
#' }
#' 
#' @export
runMethylModesApp <- function(launch.browser = TRUE, 
                              port = 3838, 
                              host = "127.0.0.1") {
  
  # Get the path to the Shiny app
  app_dir <- system.file("shiny", package = "methylModes")
  
  if (app_dir == "") {
    stop("Shiny app directory not found. Make sure the package is properly installed.")
  }
  
  # Launch the Shiny app
  shiny::runApp(app_dir, 
                launch.browser = launch.browser, 
                port = port, 
                host = host)
} 