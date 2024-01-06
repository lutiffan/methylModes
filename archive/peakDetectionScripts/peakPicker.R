library(sfsmisc)
source("/home/lutiffan/peakDetectionScripts/localMinMax.R")
# seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)

# Refined on examples:
# 173453 / "cg13332114" potentially tri-modal, with some isolated noise
# 6194 / "cg15209419" # quad-modal?
# 4283 / "cg10531774" # tall peak with an "echo" and isolated noise
# 40681 / "cg20332088" # The one that got confused between child and teen
# 185145 / "cg08274176" potentially bimodal with no gap in between
# 2167 / "cg05392448"
# 33367 / "cg19990229" # peak ratio based on height doesn't work well
# 10604 / "cg25773262" # added minDipRatio after seeing this/ why you can't just threshold at height >= 1
# 5859 / "cg14252211" # Do we skip this tiny bimodal on the left?

# 4553 Need to use peak sample volume to calculate peakPropRatio independently of peak height

# Spot checked on examples:
# 1598 / "cg04028570"
# 37109 / "cg08041448" # weird blobby trimodal

# Use for debugging: 
# source("/home/lutiffan/peakDetectionScripts/standardParams.R")
# pushToZero <- 1/ncol(betas)

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
  # seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  
  # DISTINGUISH SIGNALS FROM NOISE
  if (nrow(detected) > 1) {
    ### proportionSample check 1 ###
    leftMinx <- probeDensityEst$x[detected$leftMinIdx]
    rightMinx <- probeDensityEst$x[detected$rightMinIdx]
    
    propSample <- numeric(nrow(detected))
    for (p in seq_along(propSample)) {
      propSample[p] <- mean(betas[row.index,] > leftMinx[p] & betas[row.index,] < rightMinx[p])
    }
    detected$propSample <- propSample
    belowCutoff <- propSample < proportionSample
    ### end proportionSample check ###
    
    oneMinimaZero <- (probeDensityEst$y[detected$leftMinIdx] < pushToZero) | 
      (probeDensityEst$y[detected$rightMinIdx] < pushToZero)
    
    detected <- detected[!(belowCutoff & oneMinimaZero),]
    # seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  } else {
    detected$propSample <- 1
  }
  
  # Step 3: SPACING filter - if any peaks are too close together, compare them 
  # to their neighbors and merge shorter peaks into the tallest
  if (nrow(detected) > 1) {
    
    xValues <- probeDensityEst$x[detected$maximaIdx]
    crowded <- diff(xValues) <= personalSpace
    
    if (any(crowded)) {
      # Every other value is the number of clustered points
      clusters <- rle(crowded)$lengths + 1
      # clusteredPeakIdx contains indices for the vector detected$maximaIdx
      if (crowded[1]) {
        clusteredPeakIdx <- c(1, sort(c(which(crowded), which(diff(crowded) == 1))) + 1 )
        clusterSizes <- clusters[seq(from = 1, to = length(clusters), by = 2)]
      } else {
        clusteredPeakIdx <- sort(c(which(crowded), which(diff(crowded) == 1))) + 1; clusteredPeakIdx
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
    # seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  }
    
  if (nrow(detected) > 1) {
    ### proportionSample check 2 ###
    leftMinx <- probeDensityEst$x[detected$leftMinIdx]
    rightMinx <- probeDensityEst$x[detected$rightMinIdx]
    
    propSample <- numeric(nrow(detected))
    for (p in seq_along(propSample)) {
      propSample[p] <- mean(betas[row.index,] > leftMinx[p] & betas[row.index,] < rightMinx[p])
    }
    detected$propSample <- propSample
    belowCutoff <- propSample < proportionSample
    ### end proportionSample check ###
    detected <- detected[!belowCutoff,]
    
    # seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  }
  
  # Check the gaps or dips between peaks
  if (nrow(detected) > 1) {
    gap <- (detected$rightMinIdx[-nrow(detected)] != detected$leftMinIdx[-1]) |
      (detected$leftMinIdx[-1] < pushToZero)
    dip <- logical(nrow(detected) - 1)
    for (p in 2:nrow(detected)) {
      leftHeight <- probeDensityEst$y[detected$maximaIdx[p - 1]]
      rightHeight <- probeDensityEst$y[detected$maximaIdx[p]]
      minimumBetween <- probeDensityEst$y[detected$leftMinIdx[p]]
      # Comparing minima height to peak height
      if (leftHeight < rightHeight) {
        peakHtRatio <- probeDensityEst$y[detected$maximaIdx][p - 1]/probeDensityEst$y[detected$maximaIdx][p]
        # dipCutoff <- min(leftHeight * peakPropRatio, leftHeight * maxDipRatio)
        dipCutoff <- leftHeight * (minDipRatio + peakHtRatio * (maxDipRatio - minDipRatio))
      } else {
        peakHtRatio <- probeDensityEst$y[detected$maximaIdx][p]/probeDensityEst$y[detected$maximaIdx][p - 1]
        # dipCutoff <- min(rightHeight * peakPropRatio, rightHeight * maxDipRatio)
        dipCutoff <- rightHeight * (minDipRatio + peakHtRatio * (maxDipRatio - minDipRatio))
      }
      # # Using proportion of sample instead of height, since height approximation sometimes very imprecise
      # if (leftHeight < rightHeight) {
      #   peakPropRatio <- detected$propSample[p - 1]/detected$propSample[p]
      #   
      #   # TODO: don't tie dipCutoff to the peak volume ratio
      #   dipCutoff <- min(leftHeight * peakPropRatio, leftHeight * maxDipRatio)
      # } else {
      #   peakPropRatio <- detected$propSample[p]/detected$propSample[p - 1]
      #   dipCutoff <- min(rightHeight * peakPropRatio, rightHeight * maxDipRatio)
      # }
    
      dipCutoff <- max(dipCutoff, pushToZero)
      dip[p - 1] <- minimumBetween <= dipCutoff
    }
    
    insignificantMin <- !(gap | dip)
    
    if (any(insignificantMin)) {
      ### 
      
      # Every other value is the number of clustered points
      clusters <- rle(insignificantMin)$lengths + 1
      # clusteredPeakIdx contains indices for the vector detected$maximaIdx
      if (insignificantMin[1]) {
        clusteredPeakIdx <- c(1, sort(c(which(insignificantMin), which(diff(insignificantMin) == 1))) + 1 )
        clusterSizes <- clusters[seq(from = 1, to = length(clusters), by = 2)]
      } else {
        clusteredPeakIdx <- sort(c(which(insignificantMin), which(diff(insignificantMin) == 1))) + 1; clusteredPeakIdx
        clusterSizes <- clusters[seq(from = 2, to = length(clusters), by = 2)]
      }
      
      # Start and end indices of clusters
      clusterStart <- cumsum(clusterSizes) - clusterSizes + 1
      clusterEnd <- cumsum(clusterSizes)
      heights <- probeDensityEst$y[detected$maximaIdx[clusteredPeakIdx]]
      
      # Now we figure out which peak is the tallest in each cluster
      # Mark all insignificant peaks as "not the tallest in its cluster" by default
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
      
      # Update propSample column
      propSample <- numeric(nrow(detected))
      for (p in seq_along(propSample)) {
        propSample[p] <- mean(betas[row.index,] > leftMinx[p] & betas[row.index,] < rightMinx[p])
      }
      detected$propSample <- propSample
      
      # Sort peaks by position from left to right
      detected <- detected[order(detected$maximaIdx),]
      seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
      ###
    }
    
  }
  
  # # Try a filtering step for invariant-type distributions with noisy tails
  # if (nrow(detected) == 2) {
  #   if (probeDensityEst$y[detected$maximaIdx[1]] < probeDensityEst$y[detected$maximaIdx[2]]) {
  #     shorterIdx <- 1
  #     tallerIdx <- 2
  #   } else {
  #     shorterIdx <- 2
  #     tallerIdx <- 1
  #   }
  #   if (probeDensityEst$y[detected$maximaIdx[shorterIdx]] < 0.01 * probeDensityEst$y[detected$maximaIdx[tallerIdx]]) {
  #     detected <- detected[-shorterIdx]
  #     seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
  #   }
  # }
  
  # Step 2.5: remove small, noisy peaks area under the curve?
  # Putting this step here stops you from finding separate "echoes"
  if (length(detected$maximaIdx) > 1) { 
    densityArea <- numeric(nrow(detected))
    for (p in seq_along(densityArea)) {
      peakX <- probeDensityEst$x[detected$leftMinIdx[p]:detected$rightMinIdx[p]]
      peakY <- probeDensityEst$y[detected$leftMinIdx[p]:detected$rightMinIdx[p]]
      densityArea[p] <- integrate.xy(x = peakX, fx = peakY)
    }
    
    print(densityArea)
    #seeDetected(original.data = betas[row.index,], fitted.density = probeDensityEst, detected.peaks = detected, row.id = row.index)
    
    ## Could also remove by height, but sometimes the peaks are kind of short and wide
    # detected <- detected[probeDensityEst$y[detected$maximaIdx] > 0.5,]
  }
  
  # Step 5: check if any detected minima are close to zero, constituting a "gap"
  mins <- unique(c(detected$leftMinIdx, detected$rightMinIdx))
  if (nrow(detected) > 1) {
    gapFound <- any((detected$rightMinIdx[-nrow(detected)] != detected$leftMinIdx[-1]) |
      (detected$leftMinIdx[-1] < pushToZero))
  } else {
    gapFound <- NA
  }
  
  
  list("detected" = detected, "probeDensityEst" = probeDensityEst, "gapFound" = gapFound)
}
