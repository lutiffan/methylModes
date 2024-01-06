# USING peakSummary TABLE
peakSummaryPlot <- function(probe.names.list = NULL, row.id = NULL, probe.id = "", 
                              histogram.bins = 500, show.hypo.hyper = TRUE,
                              line.col = "orchid",
                              peak.col = "red", min.col = "blue", 
                              beta.mean.col = "magenta",
                              hypo.hyper.col = "grey") {
  noRow <- is.null(row.id)
  noProbeName <- probe.id == ""
  
  if (noRow & noProbeName) {
    missingRow <- simpleError("either row.id or probe.id is required.")
    stop(error = missingRow)
  }
  
  title <- "Beta distribution"
  # Title always include the probe name
  if (noProbeName) {
    # Retrieve probe name using row.id if probe.id not provided
    if (is.null(probe.names.list)) {
      # Assuming that betas exists in the environment
      probe.id = rownames(betas)[row.id]
    } else {
      # Not ideal to pass in a huge character vector
      probe.id = probe.names.list[row.id]
    }
    
  }
  
  title <- paste(title, "for probe", probe.id)
  
  # Row id optional. More useful for debugging
  if (!noRow) {
    title <- paste0(title, " (Row ", row.id, ")")
  } else if (is.null(probe.names.list)) {
    row.id = which(rownames(betas) == probe.id)
  } else {
    row.id = which(probe.names.list == probe.id)
  }
  
  
  # fittedDensity <- density(betas[row.id,], from = 0, to = 1, n = numBreaks, 
  #                          adjust = densityAdjust, bw = bandwidth)
  # Ignore parameter
  fittedDensity <- density(betas[row.id,], from = 0, to = 1, n = ncol(betas), 
                           adjust = densityAdjust, bw = bandwidth)
  
  # We already ran localMinMax()
  meanBeta <- peakSummary[[row.id, "meanBeta"]]
  detectedPeaks <- peakSummary[[row.id, "peakLocations"]]
  #fittedHeights <- peakSummary[[row.id, "fittedHeights"]]
  detectedMins <- c(peakSummary[[row.id, "leftMin"]], 
                    peakSummary[[row.id, "rightMin"]][peakSummary[[row.id, "numPeaks"]]])
  
  # Place the detected peaks and minima on the fitted line
  fittedPeaks <- numeric(length(detectedPeaks))
  for (p in 1:length(detectedPeaks)) {
    closestIndex <- which.min(abs(fittedDensity$x - detectedPeaks[p]))
    fittedPeaks[p] <- fittedDensity$y[closestIndex]
  }
  
  fittedMins <- numeric(length(detectedMins))
  for (m in 1:length(detectedMins)) {
    closestIndex <- which.min(abs(fittedDensity$x - detectedMins[m]))
    fittedMins[m] <- fittedDensity$y[closestIndex]
  }
  
  ggdf <- data.frame("beta" = betas[row.id,])
  
  ggplot(data = ggdf, aes(x = beta)) + 
    xlim(c(0,1)) +
    geom_histogram(bins = histogram.bins) +
    stat_density(adjust = densityAdjust, aes(color = line.col), geom = "line",
                 size = 1) +
    geom_vline(xintercept = meanBeta, aes(color = beta.mean.col)) +
    theme(legend.position = "none")
  
  hist(betas[row.id,], breaks = histogram.bins, xlim = c(0,1), 
       probability = TRUE, main = title, xlab = "Beta")
  abline(v = meanBeta, col = beta.mean.col)
  if (show.hypo.hyper) {
    abline(v = 0.3, col = hypo.hyper.col)
    abline(v = 0.7, col = hypo.hyper.col)
  }
  lines(fittedDensity, col = line.col, lwd = 2, pch = 19)
  points(detectedPeaks, fittedPeaks, 
         col = peak.col, pch = 19, cex = 1.25)
  points(detectedMins, fittedMins, 
         col = min.col, pch = 19, cex = 1.25)
  
}

peakSummaryTable <- function(probe.names.list = NULL, row.id = NULL, probe.id = "", show.boundary = T,
                             decimals = 2, var.decimals = 5) {
  noRow <- is.null(row.id)
  noProbeName <- probe.id == ""
  
  if (noRow & noProbeName) {
    missingRow <- simpleError("either row.id or probe.id is required.")
    stop(error = missingRow)
  }
  
  title <- "Beta distribution"
  # Title always include the probe name
  if (noProbeName) {
    # Retrieve probe name using row.id if probe.id not provided
    if (is.null(probe.names.list)) {
      # Assuming that betas exists in the environment
      probe.id = rownames(betas)[row.id]
    } else {
      # Not ideal to pass in a huge character vector
      probe.id = probe.names.list[row.id]
    }
    
  }
  
  title <- paste(title, "for probe", probe.id)
  
  # Row id optional. More useful for debugging
  if (!noRow) {
    title <- paste0(title, " (Row ", row.id, ")")
  } else if (is.null(probe.names.list)) {
    row.id = which(rownames(betas) == probe.id)
  } else {
    row.id = which(probe.names.list == probe.id)
  }
  
  if (show.boundary) {
    data.frame("Beta" = round(peakSummary[[row.id, "peakLocations"]], digits = decimals),
               "LeftBoundary" = round(peakSummary[[row.id, "leftMin"]], digits = decimals),
               "RightBoundary" = round(peakSummary[[row.id, "rightMin"]], digits = decimals),
               "SampleProp" = round(peakSummary[[row.id, "proportionSample"]], digits = decimals),
               "Variance" = round(peakSummary[[row.id, "peakVariance"]], digits = var.decimals))
  } else {
    data.frame("Beta" = round(peakSummary[[row.id, "peakLocations"]], digits = decimals),
               "SampleProp" = round(peakSummary[[row.id, "proportionSample"]], digits = decimals),
               "Variance" = round(peakSummary[[row.id, "peakVariance"]], digits = var.decimals))
  }
  
  
}

# GENERAL VERSION: Function to visualize detected peaks compared to histogram
# Assumes that betas exists in the environment
# Assumes that peakSummary does not exist yet
seePeaks <- function(row.id = NULL, probe.id = "", histogram.bins = 500, 
                     density.breaks = 500) {
  noRow <- is.null(row.id)
  noProbeName <- probe.id == ""
  
  if (noRow & noProbeName) {
    missingRow <- simpleError("either row.id or probe.id is required.")
    stop(error = missingRow)
  }
  
  title <- "Beta distribution"
  # Title always includes the probe name
  if (noProbeName) {
    # Retrieve probe name using row.id if probe.id not provided
    probe.id = rownames(betas)[row.id]
  }
  title <- paste(title, "for probe", probe.id)
  
  # Row id not included in title if not provided.
  if (!noRow) {
    title <- paste0(title, " (Row ", row.id, ")")
  } else {
    row.id = which(rownames(betas) == probe.id)
  }
  
  if (is.na(bandwidthType)) {
    bandwidth = "nrd0"
  } else if (bandwidthType == "sheatherJones") {
    bandwidth = bw.SJ(betas[row.id,])
  }
  fittedDensity <- density(betas[row.id,], from = 0, to = 1, n = density.breaks, 
                           adjust = densityAdjust, bw = bandwidth, 
                           kernel = kernelType)
  
  detectedPeaks <- localMinMax(fittedDensity$y, zeroThreshold = pushToZero)
  
  #### Filtering steps ####
  
  if (length(detectedPeaks$maximaIdx) > 1) {
    breaks = c(detectedPeaks$leftMinIdx, detectedPeaks$rightMinIdx[nrow(detectedPeaks)])
    #if (breaks[1] != 1) breaks <- c(1, breaks)
    #if (breaks[length(breaks)] != 500) breaks <- c(breaks, 500)
    
    detectedPeaks$propSample <- hist(betas[row.id,], breaks = fittedDensity$x[breaks], plot = FALSE)$count/length(betas[row.id,])
    aboveCutoff <- detectedPeaks$propSample > proportionSample
    detectedPeaks <- detectedPeaks[aboveCutoff,]
    
    # # DEBUG: Check plot
    # seeDetected(original.data = betas[row.id,], fitted.density = fittedDensity, detected.peaks = detectedPeaks, row.id = i)
    
  } else {
    detectedPeaks$propSample <- 1
  }
  
  # Step 4: SPACING filter - if any peaks are too close together, compare them to their neighbors and mark the shorter peak for deactivation/removal
  xValues <- fittedDensity$x[detectedPeaks$maximaIdx]
  checkCrowding <- diff(xValues) < personalSpace
  
  if (sum(checkCrowding) > 0) {
    deactivate <- rep(FALSE, nrow(detectedPeaks))
    
    for (j in 1:nrow(detectedPeaks)) {
      currentPeak <- detectedPeaks$maximaIdx[j]
      if (deactivate[j]) next # Skip if current peak already deactivated
      
      isNeighbor <- detectedPeaks$maximaIdx != currentPeak &
        xValues > xValues[j] - personalSpace & 
        xValues < xValues[j] + personalSpace
      if (sum(isNeighbor) == 0) next # No neighbors found
      
      shorter <- fittedDensity$y[detectedPeaks$maximaIdx] < fittedDensity$y[currentPeak]
      
      peaksToDeactivate <- which(isNeighbor & shorter)
      
      if (length(peaksToDeactivate) > 0) {
        suppressWarnings(deactivate[peaksToDeactivate] <- TRUE)
      }
    }
    detectedPeaks <- detectedPeaks[!deactivate,]
  }
  
  #########################
  
  hist(betas[row.id,], breaks = histogram.bins, xlim = c(0,1), probability = TRUE, main = title)
  lines(fittedDensity, col = "orchid", lwd = 2, pch = 19)
  points(fittedDensity$x[detectedPeaks$maximaIdx], fittedDensity$y[detectedPeaks$maximaIdx], col = "red", pch = 19, cex = 1.25)
  points(fittedDensity$x[detectedPeaks$leftMinIdx], fittedDensity$y[detectedPeaks$leftMinIdx], col = "blue", pch = 19, cex = 1.25)
  points(fittedDensity$x[detectedPeaks$rightMinIdx], fittedDensity$y[detectedPeaks$rightMinIdx], col = "blue", pch = 19, cex = 1.25)
}

# DEBUGGING VERSION for MAIN FOR LOOP: Function to visualize detected peaks compared to histogram
seeDetected <- function(original.data, probe.id = "", row.id = NULL, fitted.density, numBreaks = 500, detected.peaks) {
  title <- "Beta distribution"
  if (probe.id == "") {
    title <- c(title, paste("for probe", probe.id))
  }
  if (!is.null(row.id)) {
    title <- paste0("Row ", row.id)
  }
  
  hist(original.data, breaks = numBreaks, xlim = c(0,1), probability = TRUE, main = title)
  lines(fitted.density, col = "orchid", lwd = 2, pch = 19)
  points(fitted.density$x[detected.peaks$maximaIdx], fitted.density$y[detected.peaks$maximaIdx], col = "red", pch = 19, cex = 1.25)
  points(fitted.density$x[detected.peaks$leftMinIdx], fitted.density$y[detected.peaks$leftMinIdx], col = "blue", pch = 19, cex = 1.25)
  points(fitted.density$x[detected.peaks$rightMinIdx], fitted.density$y[detected.peaks$rightMinIdx], col = "blue", pch = 19, cex = 1.25)
}
