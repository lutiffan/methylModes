#' methylModes: Detection and Analysis of DNA Methylation Modes
#' 
#' @description
#' Identification of CpG sites with multimodal DNA methylation distributions, 
#' which can be caused by biological variation or technical artifacts, is an 
#' important quality control step when analyzing whole genome methylation data.
#' MethylModes is a computationally efficient, peak detection algorithm 
#' implemented in a Shiny application and easily incorporated into existing 
#' analytic pipelines for large-scale epigenetic studies.
#' 
#' @details
#' The main functions in the package are:
#' \itemize{
#'   \item \code{\link{methylModes}}: Detects methylation modes in a single probe's data
#'   \item \code{\link{methylModesBatch}}: Processes multiple probes in parallel
#'   \item \code{\link{localMinMax}}: Utility function for peak detection
#' }
#' 
#' The package also includes a Shiny application for interactive analysis,
#' accessible via \code{\link{runMethylModesApp}}.
#' 
#' @section Key Features:
#' \itemize{
#'   \item Kernel density estimation for mode detection
#'   \item Parallel processing for efficient batch analysis
#'   \item Configurable parameters for peak detection
#'   \item Interactive visualization through Shiny
#'   \item Comprehensive peak statistics and summaries
#' }
#' 
#' @section Getting Started:
#' See the vignette for a detailed tutorial:
#' \code{vignette("runMethylModesDirectly", package = "methylModes")}
#' 
#' @author
#' T. Sophia Luo \email{lutiffan@umich.edu}
#' 
#' @references
#' For more information about DNA methylation analysis:
#' \itemize{
#'   \item TODO find references
#' }
#' 
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/lutiffan/methylModes}
#'   \item Report bugs at \url{https://github.com/lutiffan/methylModes/issues}
#' }
#' 
#' @docType package
#' @name methylModes-package
NULL 

# Basic usage
runMethylModesApp()

# Custom port
runMethylModesApp(port = 8080)

# Don't open browser automatically
runMethylModesApp(launch.browser = FALSE) 