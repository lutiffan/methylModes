source("localMinMax.R")
# source("visualizeBetaPeaks.R")

# row.index = 173453 # trimodal with a short flat middle peak
# row.index = 43355 # long sparse tail

peakPicker <- function(row.index = NULL) {
  # Gets hyperparameters from the environment
  # TODO: find a smarter way to pass hyperparameters and save them all in one place
  if (is.na(bandwidthType)) {
    bandwidth = "nrd0"
  } else if (bandwidthType == "sheatherJones") {
    bandwidth = bw.SJ(betas[row.index,])
  }
  probeDensityEst <- density(betas[row.index,], from = 0, to = 1, n = numBreaks, 
                             adjust = densityAdjust, bw = bandwidth, 
                             kernel = kernelType)
  
  # Step 2: detect local maxima/minima
  detected <- localMinMax(probeDensityEst$y, zeroThreshold = pushToZero)
  #seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  
  # Step 3: SPACING filter - if any peaks are too close together, compare them 
  # to their neighbors and merge shorter peaks into the tallest
  
  zeroBefore <- probeDensityEst$y[detected$leftMinIdx[-1]] < pushToZero
  xValues <- probeDensityEst$x[detected$maximaIdx]
  tooClose <- diff(xValues) <= personalSpace
  checkCrowding <- tooClose & !zeroBefore
  
  if (sum(checkCrowding) > 0) {
    # Every other value is the number of clustered points
    clusters <- rle(checkCrowding)$lengths + 1
    # clusteredPeakIdx contains indices for the vector detected$maximaIdx
    if (checkCrowding[1]) {
      clusteredPeakIdx <- c(1, sort(c(which(checkCrowding), which(diff(checkCrowding) == 1))) + 1 )
      clusterSizes <- clusters[seq(from = 1, to = length(clusters), by = 2)]
    } else {
      clusteredPeakIdx <- sort(c(which(checkCrowding), which(diff(checkCrowding) == 1))) + 1; clusteredPeakIdx
      clusterSizes <- clusters[seq(from = 2, to = length(clusters), by = 2)]
    }
    
    # Start and end indices (for the vector clusteredPeakIdx) clusters
    clusterStart <- cumsum(clusterSizes) - clusterSizes + 1
    clusterEnd <- cumsum(clusterSizes)
    heights <- probeDensityEst$y[detected$maximaIdx[clusteredPeakIdx]]
    
    # Now we figure out which peak is the tallest in each cluster
    # Mark all crowded peaks as "not the tallest in its cluster" by default
    clusterTallest <- logical(length(clusteredPeakIdx))
    # Get the tallest peak in each cluster
    for (x in seq_along(clusterSizes)) {
      clusterWinner <- which.max(heights[clusterStart[x]:clusterEnd[x]]) + 
        clusterStart[x] - 1
      
      clusterTallest[clusterWinner] <- TRUE
      detected$leftMinIdx[clusteredPeakIdx[clusterWinner]] <- detected$leftMinIdx[clusteredPeakIdx[clusterStart[x]]]
      detected$rightMinIdx[clusteredPeakIdx[clusterWinner]] <- detected$rightMinIdx[clusteredPeakIdx[clusterEnd[x]]]
    }
    
    clearWinners <- setdiff(1:nrow(detected), clusteredPeakIdx)
    
    detected <- detected[c(clearWinners, clusteredPeakIdx[clusterTallest]),]
    
    # Sort peaks by position from left to right
    detected <- detected[order(detected$maximaIdx),]
  }
  # seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  
  # Step 4: SIZE filter - (only if multiple maxima found) remove peaks that represent too small of a proportion of the sample to be significant
  # if (length(detected$maximaIdx) > 1) {
  #   breaks <- c(detected$leftMinIdx, detected$rightMinIdx[nrow(detected)])
  # 
  #   detected$propSample <- hist(betas[row.index,], breaks = probeDensityEst$x[breaks], plot = FALSE)$count/length(betas[row.index,])
  #   aboveCutoff <- detected$propSample > proportionSample
  #   detected <- detected[aboveCutoff,]
  # 
  # } else {
  #   detected$propSample <- 1
  # }
  
  if (length(detected$maximaIdx) > 1) {
    breaks <- c(detected$leftMinIdx, detected$rightMinIdx[nrow(detected)])
    
    detected$propSample <- hist(betas[row.index,], breaks = probeDensityEst$x[breaks], plot = FALSE)$count/length(betas[row.index,])
    belowCutoff <- detected$propSample < proportionSample
    peaksToMerge <- which(belowCutoff)
    for (b in seq_along(peaksToMerge)) {
      currentPeak <- peaksToMerge[b]
      if (detected$leftMinIdx[b] > detected$rightMinIdx[b]) { # Merge with left neighbor
        
      } else if (detected$leftMinIdx[b] < detected$rightMinIdx[b]) { # Merge with right neighbor
        
      } else { # Tiebreaker 1
        sum(detected$propSample[1:b]) 
        
      }
    }
    
    detected <- detected[!belowCutoff,]
    
  } else {
    detected$propSample <- 1
  }
  
  # Step 5: check if any detected minima are close to zero, constituting a "gap"
  mins <- c(detected$leftMinIdx, detected$rightMinIdx[nrow(detected)])
  if (nrow(detected) > 1) {
    gap <- sum(probeDensityEst$y[mins[2:(length(mins) - 1)]] <= pushToZero) > 0
  } else {
    gap <- FALSE
  }
  
  
  list("detected" = detected, "probeDensityEst" = probeDensityEst, "gap" = gap)
}
