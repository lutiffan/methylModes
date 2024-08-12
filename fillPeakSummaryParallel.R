# Requires that hyperparameters exist in the environment
fillPeakSummaryParallel <- function(betas = NULL) {

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
  
  # .export = ls(envir = globalenv())
  
  
  peakSummary <- foreach::foreach(probe = iter(betas, by = "row"), 
                                  .combine = "rbind",
                                  .packages = c("foreach", "data.table"),
                                  .export = c("methylModes",
                                              "localMinMax",
                                              "proportionSample",
                                              "peakDistance",
                                              "kernelType",
                                              "bandwidthType",
                                              "numBreaks",
                                              "densityAdjust",
                                              "pushToZero")) %dopar% {
    
    # Steps 1-4: smooth histogram, detect local maxima/minima, filter by spacing, filter by sample %, detect presence of any "gaps"
    foundPeaks <- methylModes(row.data = probe)
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
  stopCluster(cl)
  
  peakSummary
}
