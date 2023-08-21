source("peakPicker.R")

# Requires that peakSummary and associated hyperparameters exist in the environment
fillPeakSummary <- function() {
  for (i in 1:totalRows) { #nrow(betas)) {
    rowIndex <- rangeStart + i - 1
    # Step 1: smooth the histogram
    # Setting bandwidth parameter to bw.SJ() seems to fit data very well
    # but may be inefficient and it is noisy
    
    # Steps 2-4: detect local minima, filter by spacing, filter by sample %
    foundPeaks <- peakPicker(row.index = rowIndex)
    detected <- foundPeaks$detected
    probeDensityEst <- foundPeaks$probeDensityEst
    
    # Step 5: Calculate summary statistics
    # 5a: variance for observations included in each peak
    peakVariance = numeric(nrow(detected))
    for (k in 1:nrow(detected)) {
      leftBeta <- probeDensityEst$x[detected$leftMinIdx[k]]
      rightBeta <- probeDensityEst$x[detected$rightMinIdx[k]]
      peakVariance[k] <- var(betas[rowIndex,][betas[rowIndex,] > leftBeta & betas[rowIndex,] < rightBeta])
    }
    
    # 5b: mean beta value among all samples
    meanBeta <- mean(betas[rowIndex,])
    
    peakSummary[i,] <- list("numPeaks" = nrow(detected), 
                            "meanBeta" = meanBeta,
                            "peakLocations" = probeDensityEst$x[detected$maximaIdx],
                            # "fittedHeights" = probeDensityEst$y[detected$maximaIdx],
                            "leftMin" = probeDensityEst$x[detected$leftMinIdx],
                            "rightMin" = probeDensityEst$x[detected$rightMinIdx],
                            "proportionSample" = detected$propSample,
                            "peakVariance" = peakVariance)
  }
  peakSummary
}

# iterateOverProbes <- function() {
#   # Iterate over probes
#   startTime <- Sys.time()
#   for (i in 1:totalRows) { #nrow(betas)) {
#     rowIndex <- rangeStart + i - 1
#     # Step 1: smooth the histogram
#     # Setting bandwidth parameter to bw.SJ() seems to fit data well
#     # but may not be efficient
#     if (is.na(bandwidthType)) {
#       bandwidth = "nrd0"
#     } else if (bandwidthType == "sheatherJones") {
#       bandwidth = bw.SJ(betas[rowIndex,])
#     }
#     probeDensityEst <- density(betas[rowIndex,], from = 0, to = 1, n = numBreaks, 
#                                adjust = densityAdjust, bw = bandwidth, 
#                                kernel = kernelType)
#     
#     # Step 2: detect local maxima/minima
#     detected <- localMinMax(probeDensityEst$y, zeroThreshold = pushToZero)
#     # Check plot for debugging
#     # seeDetected(original.data = betas[rowIndex,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = i)
#     
#     # Step 3: SIZE filter - (only if multiple maxima found) remove peaks that represent too small of a proportion of the sample to be significant
#     if (length(detected$maximaIdx) > 1) {
#       breaks = c(detected$leftMinIdx, detected$rightMinIdx[nrow(detected)])
#       #if (breaks[1] != 1) breaks <- c(1, breaks)
#       #if (breaks[length(breaks)] != 500) breaks <- c(breaks, 500)
#       
#       detected$propSample <- hist(betas[rowIndex,], breaks = probeDensityEst$x[breaks], plot = FALSE)$count/length(betas[rowIndex,])
#       aboveCutoff <- detected$propSample > proportionSample
#       detected <- detected[aboveCutoff,]
#       
#       # # DEBUG: Check plot
#       # seeDetected(original.data = betas[rowIndex,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = i)
#       
#     } else {
#       detected$propSample <- 1
#     }
#     
#     # Step 4: SPACING filter - if any peaks are too close together, compare them to their neighbors and mark the shorter peak for deactivation/removal
#     xValues <- probeDensityEst$x[detected$maximaIdx]
#     checkCrowding <- diff(xValues) < personalSpace
#     
#     if (sum(checkCrowding) > 0) {
#       deactivate <- rep(FALSE, nrow(detected))
#       
#       for (j in 1:nrow(detected)) {
#         currentPeak <- detected$maximaIdx[j]
#         if (deactivate[j]) next # Skip if current peak already deactivated
#         
#         isNeighbor <- detected$maximaIdx != currentPeak &
#           xValues > xValues[j] - personalSpace & 
#           xValues < xValues[j] + personalSpace
#         if (sum(isNeighbor) == 0) next # No neighbors found
#         
#         shorter <- probeDensityEst$y[detected$maximaIdx] < probeDensityEst$y[currentPeak]
#         
#         peaksToDeactivate <- which(isNeighbor & shorter)
#         
#         if (length(peaksToDeactivate) > 0) {
#           suppressWarnings(deactivate[peaksToDeactivate] <- TRUE)
#         }
#       }
#       detected <- detected[!deactivate,]
#     }
#     
#     # # DEBUG: Check plot
#     # seeDetected(original.data = betas[rowIndex,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = i)
#     
#     # # More debug code
#     # print(detected)
#     # browser()
#     
#     # Step 5: Calculate summary statistics
#     # 5a: variance for observations included in each peak
#     peakVariance = numeric(nrow(detected))
#     for (k in 1:nrow(detected)) {
#       leftBeta <- probeDensityEst$x[detected$leftMinIdx[k]]
#       rightBeta <- probeDensityEst$x[detected$rightMinIdx[k]]
#       peakVariance[k] <- var(betas[rowIndex,][betas[rowIndex,] > leftBeta & betas[rowIndex,] < rightBeta])
#     }
#     
#     # 5b: invariance (mean beta is less than 0.3 or greater than 0.7) 
#     meanBeta <- mean(betas[rowIndex,])
#     
#     peakSummary[i,] <- list("numPeaks" = nrow(detected), 
#                             "meanBeta" = meanBeta,
#                             "peakLocations" = probeDensityEst$x[detected$maximaIdx],
#                             # "fittedHeights" = probeDensityEst$y[detected$maximaIdx],
#                             "leftMin" = probeDensityEst$x[detected$leftMinIdx],
#                             "rightMin" = probeDensityEst$x[detected$rightMinIdx],
#                             "proportionSample" = detected$propSample,
#                             "peakVariance" = peakVariance)
#   }
#   endTime <- Sys.time() - startTime
#   
#   sink("sensitivityAnalysis/runtimes.txt", append = TRUE)
#   print(output)
#   print(endTime)
#   sink()
#   
#   saveRDS(peakSummary, file = output)
# }
