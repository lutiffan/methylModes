source("localMinMax.R")

# Clean slate for testing : rm(list = setdiff(ls(), "betas"))
# debug example: row.index = 4661
# Clusters: 1 2 3 4, 5 6 7, 10 11
# Winners: 1, 7, 8

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
  # Check plot for debugging
  # seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  
  
  # Step 3: SPACING filter - if any peaks are too close together, compare them 
  # to their neighbors and merge shorter peaks into the tallest
  xValues <- probeDensityEst$x[detected$maximaIdx]
  checkCrowding <- diff(xValues) <= personalSpace
  
  if (sum(checkCrowding) > 0) {
    # Count how many peaks are in each cluster
    clusters <- rle(checkCrowding)$lengths + 1
    # clusters is a vector of indices for the vector detected$maximaIdx
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
    
    #### Wrong ####
    #heights <- probeDensityEst$y[detected$maximaIdx]
    #### Wrong ####
    
    # Keep track of indices (of maximaIdx values) marking start/end of the true peaks
    ### Wrong ###
    #numClusters <- sum(diff(clusteredPeakIdx) >= 2) + 1
    ###       ###
    
    
    # numClusters <- length(clusterSizes)
    # bookEnds <- data.frame("start" = integer(numClusters), "end" = integer(numClusters))
    
    # # Initialize using first crowded peak
    # bookEnds$start[1] <- clusteredPeakIdx[1]
    
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
    
    # if (numClusters == 1) {
    #   clusterTallest[which.max(heights[clusteredPeakIdx])] <- TRUE
    #   # TODO: take start and end minima of cluster
    #   # bookEnds$end[1] <- clusteredPeakIdx[length(clusteredPeakIdx)]
    # } else {
    #   # Initialize for loop variables using first crowded peak
    #   clusterIdx <- 1
    #   clusterTallestIdx <- 1
    #   clusterTallestHt <- heights[clusteredPeakIdx[1]]; clusterTallestHt
    #   numCrowded <- length(clusteredPeakIdx)
    #   
    #   #### Wrong ####
    #   for (x in 2:numCrowded) {
    #     if ((clusteredPeakIdx[x] - clusteredPeakIdx[x - 1]) > 1) {
    #       # debug
    #       # print("new cluster")
    #       # Mark the end of the old cluster
    #       bookEnds$end[clusterIdx] <- clusteredPeakIdx[x - 1]
    # 
    #       # Mark the previous cluster's tallest
    #       clusterTallest[clusterTallestIdx] <- TRUE
    # 
    #       clusterTallestIdx <- x
    #       clusterTallestHt <- heights[clusteredPeakIdx[x]]
    # 
    #       clusterIdx <- clusterIdx + 1
    # 
    #       # Keep track of the start of the new cluster
    #       bookEnds$start[clusterIdx] <- clusteredPeakIdx[x]
    #     }
    # 
    #     if (heights[clusteredPeakIdx[x]] > clusterTallestHt) {
    #       # debug
    #       # print("new tallest peak")
    #       clusterTallestIdx <- x
    #       clusterTallestHt <- heights[clusteredPeakIdx[x]]
    #     }
    # 
    #     if (x == numCrowded) {
    #       # Mark the end of the last cluster
    #       bookEnds$end[clusterIdx] <- clusteredPeakIdx[numCrowded]
    # 
    #       clusterTallest[clusterTallestIdx] <- TRUE
    #     }
    #     # debug - remove before running whole for loop!
    #     # x = x + 1; x
    #     # clusterTallestIdx; clusterTallestHt
    #   }
    #   #### Wrong ####
    # }
    
    # newPeaks <- clusteredPeakIdx[clusterTallest]
    
    # for (x in seq_along(newPeaks)) {
    #   detected$leftMinIdx[newPeaks[x]] <- detected$leftMinIdx[bookEnds$start[x]]
    #   detected$rightMinIdx[newPeaks[x]] <- detected$rightMinIdx[bookEnds$end[x]]
    # }
    
    clearWinners <- setdiff(1:nrow(detected), clusteredPeakIdx)
    
    detected <- detected[c(clearWinners, clusteredPeakIdx[clusterTallest]),]
    
    # Sort peaks by position from left to right
    detected <- detected[order(detected$maximaIdx),]
  }
  # seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  
  # Step 4: SIZE filter - (only if multiple maxima found) remove peaks that represent too small of a proportion of the sample to be significant
  if (length(detected$maximaIdx) > 1) {
    breaks = c(detected$leftMinIdx, detected$rightMinIdx[nrow(detected)])
    #if (breaks[1] != 1) breaks <- c(1, breaks)
    #if (breaks[length(breaks)] != 500) breaks <- c(breaks, 500)
    
    detected$propSample <- hist(betas[row.index,], breaks = probeDensityEst$x[breaks], plot = FALSE)$count/length(betas[row.index,])
    aboveCutoff <- detected$propSample > proportionSample
    detected <- detected[aboveCutoff,]
    
    # # DEBUG: Check plot
    # seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
    
  } else {
    detected$propSample <- 1
  }
  
  # # DEBUG: Check plot
  # print(seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = index))
  
  # # More debug code
  # print(detected)
  #browser()
  list("detected" = detected, "probeDensityEst" = probeDensityEst)
}
