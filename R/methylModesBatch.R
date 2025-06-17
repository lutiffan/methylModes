#' Batch Process Methylation Mode Detection
#' 
#' This function processes multiple probes (rows) of DNA methylation data in 
#' parallel to detect modes in methylation data for each probe. It uses parallel 
#' processing to efficiently handle large datasets and returns a summary of the 
#' detected peaks for each probe.
#' 
#' @param betas A matrix (do not use a data.frame) of methylation beta values 
#'   where rows are probes and columns are samples. Row names should be probe IDs.
#' @param proportionSample Numeric threshold for minimum proportion of samples 
#'   in a peak (default: 0.01).
#' @param peakDistance Numeric minimum distance between peaks 
#'   (default: 0.10).
#' @param kernelType Character string specifying the kernel type 
#'   (default: "gaussian").
#' @param bandwidthType Character string specifying the bandwidth selection 
#'   method (default: NA).
#' @param numBreaks Integer number of points for density estimation 
#'   (default: 500).
#' @param densityAdjust Numeric adjustment factor for density estimation 
#'   (default: 1.5).
#' @param pushToZero Numeric threshold for pushing small density values to zero
#'   (default: 1e-6).
#' 
#' @return A data.table with one row per probe containing:
#' \describe{
#'   \item{probeName}{Character vector of probe IDs}
#'   \item{numPeaks}{Numeric vector of number of peaks detected per probe}
#'   \item{meanBeta}{Numeric vector of mean beta value per probe}
#'   \item{peakLocations}{List of numeric vectors containing peak locations 
#'         for each probe}
#'   \item{leftMin}{List of numeric vectors containing left minima for each peak}
#'   \item{rightMin}{List of numeric vectors containing right minima for each peak}
#'   \item{proportionSample}{List of numeric vectors containing proportion of 
#'         samples in each peak}
#'   \item{peakVariance}{List of numeric vectors containing variance of values 
#'         in each peak}
#' }
#' 
#' @details
#' The function processes the data in parallel using the following steps:
#' \enumerate{
#'   \item Creates a parallel cluster using available cores (minus one)
#'   \item Processes each probe in parallel using methylModes()
#'   \item Combines results into a data.table
#'   \item Calculates additional statistics (mean, variance) for each peak
#' }
#' 
#' @examples
#' # Generate example methylation data
#' set.seed(123)
#' probes <- 10
#' samples <- 100
#' betas <- matrix(runif(probes * samples), nrow = probes, ncol = samples)
#' rownames(betas) <- paste0("cg", 1:probes)
#' 
#' # Run batch processing with default parameters
#' results <- methylModesBatch(betas)
#' 
#' # Run with custom parameters
#' custom_params <- list(proportionSample = 0.1, peakDistance = 0.2)
#' results <- methylModesBatch(betas, proportionSample = 0.1, peakDistance = 0.2)
#' 
#' # View summary of results
#' print(results[, .(probeName, numPeaks, meanBeta)])
#' 
#' # Plot distribution of peak counts
#' hist(results$numPeaks, main = "Distribution of Peak Counts",
#'      xlab = "Number of Peaks", ylab = "Frequency")
#' 
#' @export
methylModesBatch <- function(betas = NULL,
                            proportionSample = 0.01,
                            peakDistance = 0.10,
                            kernelType = "gaussian",
                            bandwidthType = NA,
                            numBreaks = 500,
                            densityAdjust = 1.5,
                            pushToZero = 1e-6) {

  if (is.null(betas)) stop(simpleError("Invalid beta matrix."))
  
  template <- data.table("probeName" = character(1),
                         "numPeaks" = numeric(1),
                         "meanBeta" = numeric(1),
                         "peakLocations" = vector(mode = "list", 
                                                  length = 1),
                         "leftMin" = vector(mode = "list", 
                                            length = 1),
                         "rightMin" = vector(mode = "list", 
                                             length = 1),
                         "proportionSample" = vector(mode = "list", 
                                                     length = 1),
                         "peakVariance" = vector(mode = "list", 
                                                 length = 1))
  
  # Make new cluster
  cl <- parallelly::makeClusterPSOCK(availableCores(omit = 1), autoStop = TRUE)
  
  # Ensure the cluster is stopped even if the function exits because of an error
  # on.exit(stopCluster(cl), add = TRUE)
  
  # Register parallel backend
  doParallel::registerDoParallel(cl)
  
  peakSummary <- foreach::foreach(probe = iter(betas, by = "row"), 
                                  .combine = "rbind",
                                  .packages = c("foreach", "data.table", "methylModes")) %dopar% {
    
    # Steps 1-4: smooth histogram, detect local maxima/minima, filter by 
    # spacing, filter by sample %, detect presence of any "gaps"
    foundPeaks <- methylModes(row.data = probe,
                              proportionSample = proportionSample,
                              peakDistance = peakDistance,
                              kernelType = kernelType,
                              bandwidthType = bandwidthType,
                              numBreaks = numBreaks,
                              densityAdjust = densityAdjust,
                              pushToZero = pushToZero)
    detected <- foundPeaks$detected
    probeDensityEst <- foundPeaks$probeDensityEst
    
    # Step 5: Calculate summary statistics
    # 5a: variance for observations included in each peak
    peakVariance = numeric(nrow(detected))
    for (k in 1:nrow(detected)) {
      leftBeta <- probeDensityEst$x[detected$leftMinIdx[k]]
      rightBeta <- probeDensityEst$x[detected$rightMinIdx[k]]
      peakVariance[k] <- var(probe[((probe > leftBeta) & (probe < rightBeta))])
    }
    
    # 5b: mean beta value among all samples
    meanBeta <- mean(probe)
    
    template[1,] <- list("probeName" = rownames(probe),
      "numPeaks" = nrow(detected), 
      "meanBeta" = meanBeta,
      "peakLocations" = probeDensityEst$x[detected$maximaIdx],
      # "fittedHeights" = probeDensityEst$y[detected$maximaIdx],
      "leftMin" = probeDensityEst$x[detected$leftMinIdx],
      "rightMin" = probeDensityEst$x[detected$rightMinIdx],
      "proportionSample" = detected$propSample,
      "peakVariance" = peakVariance) #,
      # "gapFound" = foundPeaks$gap) 
    # The complicated version of methylModes has foundPeaks$gapFound
    template
  }
  # Close connections (not strictly necessary, but best practice)
  parallel::stopCluster(cl)
  
  peakSummary
}
