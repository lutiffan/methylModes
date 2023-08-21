library(sfsmisc)
source("/home/lutiffan/peakDetectionScripts/localMinMax.R")
# seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)

# Refined on examples:
# 173453 / "cg13332114"
# 6194 / "cg15209419"
# 4283 / "cg10531774"

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
  
  # # Step 2.5: SIZE filter - (only if multiple maxima found) remove peaks that represent too small of a proportion of the sample to be significant
  # if (length(detected$maximaIdx) > 1) {
  #   isolated <- (probeDensityEst$y[detected$leftMinIdx] < pushToZero) & (probeDensityEst$y[detected$rightMinIdx] < pushToZero)
  #   breaks = c(detected$leftMinIdx, detected$rightMinIdx[nrow(detected)])
  #   
  #   detected$propSample <- hist(betas[row.index,], breaks = probeDensityEst$x[breaks], plot = FALSE)$count/length(betas[row.index,])
  #   noiseIsland <- (detected$propSample < proportionSample) & isolated
  #   detected <- detected[!noiseIsland,]
  #   
  # } else {
  #   detected$propSample <- 1
  # }
  
  #seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  
  # Step 2.5: remove small, noisy peaks area under the curve?
  # Putting this step here stops you from finding separate "echoes"
  if (length(detected$maximaIdx) > 1) { 
    densityArea <- numeric(nrow(detected))
    for (p in seq_along(densityArea)) {
      peakX <- probeDensityEst$x[detected$leftMinIdx[p]:detected$rightMinIdx[p]]
      peakY <- probeDensityEst$y[detected$leftMinIdx[p]:detected$rightMinIdx[p]]
      densityArea[p] <- integrate.xy(x = peakX, fx = peakY)
    }
    
    detected <- detected[densityArea > 0.05,]
    #seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
    
    ## Could also remove by height, but sometimes the peaks are kind of short and wide
    # detected <- detected[probeDensityEst$y[detected$maximaIdx] > 0.5,]
  }
  
  # Step 3: SPACING filter - if any peaks are too close together, compare them 
  # to their neighbors and merge shorter peaks into the tallest
  xValues <- probeDensityEst$x[detected$maximaIdx]
  checkCrowding <- diff(xValues) <= personalSpace
  
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
    
    # Start and end indices of clusters
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
  #seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  
  # Step 4: SIZE filter - (only if multiple maxima found) remove peaks that represent too small of a proportion of the sample to be significant
  if (length(detected$maximaIdx) > 1) { 
    leftMinx <- probeDensityEst$x[detected$leftMinIdx]
    rightMinx <- probeDensityEst$x[detected$rightMinIdx]
    
    propSample <- numeric(nrow(detected))
    for (p in seq_along(propSample)) {
      propSample[p] <- mean(betas[row.index,] > leftMinx[p] & betas[row.index,] < rightMinx[p])
    }
    detected$propSample <- propSample
    aboveCutoff <- propSample > proportionSample
    detected <- detected[aboveCutoff,]
  } else {
      detected$propSample <- 1
  }
  
  # if (length(detected$maximaIdx) > 1) {
  #   mins <- unique(c(detected$leftMinIdx, detected$rightMinIdx))
  #   breaks <- probeDensityEst$x[mins]
  #   #breaks = c(detected$leftMinIdx, detected$rightMinIdx[nrow(detected)])
  #   peakCounts <- hist(betas[row.index,], breaks = breaks, plot = FALSE)$count
  #   # wrong # detected$propSample <- peakCounts[seq(1, length(peakCounts), by = 2)]/length(betas[row.index,])
  #   aboveCutoff <- detected$propSample > proportionSample
  #   detected <- detected[aboveCutoff,]
  # }
  #seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  
  # Step 5: check if any detected minima are close to zero, constituting a "gap"
  mins <- unique(c(detected$leftMinIdx, detected$rightMinIdx))
  if (nrow(detected) > 1) {
    gap <- sum(probeDensityEst$y[mins[2:(length(mins) - 1)]] < pushToZero) > 0
  } else {
    gap <- NA
  }
  
  
  list("detected" = detected, "probeDensityEst" = probeDensityEst, "gap" = gap)
}
