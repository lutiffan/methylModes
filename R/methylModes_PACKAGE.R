#' methylModes: Detection of Multi-Modal Distributions in DNA Methylation Data
#' 
#' @description
#' Identification of CpG sites with multimodal DNA methylation distributions, 
#' which can be caused by biological variation or technical artifacts, is an 
#' important quality control step when analyzing whole genome methylation data. 
#' MethylModes is a computationally efficient, peak detection algorithm implemented 
#' in a Shiny application and easily incorporated into existing analytic pipelines 
#' for large-scale epigenetic studies.
#' 
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{methylModes}}: Detect methylation modes in a single probe
#'   \item \code{\link{methylModesBatch}}: Batch process multiple probes
#'   \item \code{\link{runMethylModesApp}}: Launch the Shiny application
#'   \item \code{\link{localMinMax}}: Detect local minima and maxima
#' }
#' 
#' @section Datasets:
#' \itemize{
#'   \item \code{\link{IlluminaManifest450k}}: 450k array annotation data
#'   \item \code{\link{IlluminaManifestEPIC}}: EPIC array annotation data
#'   \item \code{\link{toyMultimodalData}}: Test dataset for examples
#' }
#' 
#' @section Shiny Application:
#' The package includes a comprehensive Shiny application for interactive analysis:
#' \itemize{
#'   \item Upload and validate methylation data
#'   \item Interactive visualization of results
#'   \item Batch processing capabilities
#'   \item Export results in various formats
#' }
#' 
#' @author T. Sophia Luo \email{lutiffan@umich.edu}
#' @author Jonathon LeFaive
#' @author John Dou
#' @author Kelly M. Bakulski
#' @author Erin B. Ware
#' @author Matthew Zawistowski
#' 
#' @references
#' For more information, see the package vignette: \code{vignette("methylModes")}
#' 
#' @keywords internal
"_PACKAGE"
