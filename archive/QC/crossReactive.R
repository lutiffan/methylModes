# !!! First I'm running the first two chunks of forLoop_peakDetection.Rmd

load("/nfs/turbo/bakulski1/cross.probes.info.450k.rda")
head(cross.probes.info)

# Get vector of names of cross-reactive probes in our beta matrix
crRows <- which(rownames(betas) %in% cross.probes.info$TargetID)

rangeStart = 1
rangeEnd = length(crRows)
totalRows = rangeEnd - rangeStart + 1

readRDS("crossReactivePeakSummary.RDS")

ggplot(data = peakSummary, aes(numPeaks)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust = -1) +
  ggtitle(label = paste0("Cross-reactive probe summary ", 
                         "(Prop. = ", proportionSample, 
                         " Space = ", personalSpace, ")"))

##### Run peak detection on cross-reactive probes #####

# Container for all probes
peakSummary <- data.table("numPeaks" = numeric(totalRows), 
                          "meanBeta" = numeric(totalRows),
                          "peakLocations" = vector(mode = "list", 
                                                   length = totalRows),
                          # "fittedHeights" = vector(mode = "list", 
                          #                     length = totalRows),
                          "leftMin" = vector(mode = "list", 
                                             length = totalRows),
                          "rightMin" = vector(mode = "list", 
                                              length = totalRows),
                          "proportionSample" = vector(mode = "list", 
                                                      length = totalRows),
                          "peakVariance" = vector(mode = "list", 
                                                  length = totalRows))

startTime <- Sys.time()
for (i in 1:totalRows) { #nrow(betas)) {
  rowIndex <- crRows[i]
  # Step 1: smooth the histogram
  # Setting bandwidth parameter to bw.SJ() seems to fit data well
  # but may not be efficient
  if (is.na(bandwidthType)) {
    bandwidth = "nrd0"
  } else if (bandwidthType == "sheatherJones") {
    bandwidth = bw.SJ(betas[rowIndex,])
  }
  probeDensityEst <- density(betas[rowIndex,], from = 0, to = 1, n = numBreaks, 
                             adjust = densityAdjust, bw = bandwidth, 
                             kernel = kernelType)
  
  # Step 2: detect local maxima/minima
  detected <- localMinMax(probeDensityEst$y, zeroThreshold = pushToZero)
  # Check plot for debugging
  # seeDetected(original.data = betas[rowIndex,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = i)
  
  # Step 3: SIZE filter - (only if multiple maxima found) remove peaks that represent too small of a proportion of the sample to be significant
  if (length(detected$maximaIdx) > 1) {
    breaks = c(detected$leftMinIdx, detected$rightMinIdx[nrow(detected)])
    #if (breaks[1] != 1) breaks <- c(1, breaks)
    #if (breaks[length(breaks)] != 500) breaks <- c(breaks, 500)
    
    detected$propSample <- hist(betas[rowIndex,], breaks = probeDensityEst$x[breaks], plot = FALSE)$count/length(betas[rowIndex,])
    aboveCutoff <- detected$propSample > proportionSample
    detected <- detected[aboveCutoff,]
    
    # # DEBUG: Check plot
    # seeDetected(original.data = betas[rowIndex,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = i)
    
  } else {
    detected$propSample <- 1
  }
  
  # Step 4: SPACING filter - if any peaks are too close together, compare them to their neighbors and mark the shorter peak for deactivation/removal
  xValues <- probeDensityEst$x[detected$maximaIdx]
  checkCrowding <- diff(xValues) < personalSpace
  
  if (sum(checkCrowding) > 0) {
    deactivate <- rep(FALSE, nrow(detected))
    
    for (j in 1:nrow(detected)) {
      currentPeak <- detected$maximaIdx[j]
      if (deactivate[j]) next # Skip if current peak already deactivated
      
      isNeighbor <- detected$maximaIdx != currentPeak &
        xValues > xValues[j] - personalSpace & 
        xValues < xValues[j] + personalSpace
      if (sum(isNeighbor) == 0) next # No neighbors found
      
      shorter <- probeDensityEst$y[detected$maximaIdx] < probeDensityEst$y[currentPeak]
      
      peaksToDeactivate <- which(isNeighbor & shorter)
      
      if (length(peaksToDeactivate) > 0) {
        suppressWarnings(deactivate[peaksToDeactivate] <- TRUE)
      }
    }
    detected <- detected[!deactivate,]
  }
  
  # # DEBUG: Check plot
  # seeDetected(original.data = betas[rowIndex,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = i)
  
  # # More debug code
  # print(detected)
  # browser()
  
  # Step 5: Calculate summary statistics
  # 5a: variance for observations included in each peak
  peakVariance = numeric(nrow(detected))
  for (k in 1:nrow(detected)) {
    leftBeta <- probeDensityEst$x[detected$leftMinIdx[k]]
    rightBeta <- probeDensityEst$x[detected$rightMinIdx[k]]
    peakVariance[k] <- var(betas[rowIndex,][betas[rowIndex,] > leftBeta & betas[rowIndex,] < rightBeta])
  }
  
  # 5b: invariance (mean beta is less than 0.3 or greater than 0.7) 
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
endTime <- Sys.time() - startTime; endTime

# saveRDS(peakSummary, "crossReactivePeakSummary_like_gaphunter.RDS")
