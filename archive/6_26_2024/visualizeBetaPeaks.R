source("/home/lutiffan/peakDetectionScripts/methylModes.R")

# USING peakSummary AND BETA MATRIX
peakSummaryPlot <- function(beta.data, peak.summary, probe.id = NULL, 
                              # range.start = NULL, range.end = NULL,
                              histogram.bins = 50, histogram.col = "black",
                              show.hypo.hyper = TRUE,
                              line.col = "orchid",
                              peak.col = "red", min.col = "blue", 
                              beta.mean.col = "magenta",
                              hypo.hyper.col = "grey",
                              label.col = "red") {
  
  # Try to extract missing probe.id from provided data 
  # (assuming it was subsetted from a matrix with drop = FALSE)
  if (is.null(probe.id)) {
    probe.id <- rownames(beta.data)
  }
  
  # If that failed, we can't label the plot with the probe name
  if (is.null(probe.id)) {
    missingRow <- simpleError("ensure that you set the 'drop' argument to FALSE when passing in a subsetted matrix row, or that you provide probe.id.")
    stop(error = missingRow)
  }
  
  title <- paste("Beta distribution for probe", probe.id)
  
  # # Row id optional in title.
  # if (!noRow) { # probe.id must have been provided
  #   title <- paste0(title, " (peakSummary row ", row.id, ")")
  # } else {
  #   row.id = which(rownames(betas) == probe.id)
  # } 
  
  # betaIndex <- range.start + row.id - 1
  
  if (is.na(bandwidthType)) {
    bandwidth = "nrd0"
  } else if (bandwidthType == "sheatherJones") {
    bandwidth = bw.SJ(betas[row.id,])
  }
  
  fittedDensity <- density(beta.data, from = 0, to = 1, n = numBreaks, 
                           adjust = densityAdjust, bw = bandwidth)
  
  # We already ran localMinMax()
  meanBeta <- peak.summary$meanBeta
  detectedPeaks <- unlist(peak.summary$peakLocations)
  #fittedHeights <- peak.summary[[, "fittedHeights"]]
  detectedMins <- c(unlist(peak.summary$leftMin), 
                    unlist(peak.summary$rightMin[peak.summary$numPeaks]))
  
  # Place the detected peaks and minima on the fitted line
  fittedPeaks <- numeric(length(detectedPeaks))
  for (p in 1:length(detectedPeaks)) {
    closestIndex <- which.min(abs(fittedDensity$x - detectedPeaks[p]))
    fittedPeaks[p] <- fittedDensity$y[closestIndex]
  }
  print("Fitted peak heights:")
  print(fittedPeaks)
  
  fittedMins <- numeric(length(detectedMins)) 
  for (m in 1:length(detectedMins)) {
    closestIndex <- which.min(abs(fittedDensity$x - detectedMins[m]))
    fittedMins[m] <- fittedDensity$y[closestIndex]
  }
  
  print("Fitted min heights:")
  print(fittedMins)
  
  print("Proportion of sample in each peak, from left to right:")
  print(peak.summary$proportionSample)
  
  hist(beta.data, breaks = histogram.bins, xlim = c(0,1), 
       probability = TRUE, main = title, xlab = "Beta", col = histogram.col)
  abline(v = meanBeta, col = beta.mean.col)
  if (show.hypo.hyper) {
    abline(v = 0.3, col = hypo.hyper.col)
    abline(v = 0.7, col = hypo.hyper.col)
  }
  lines(fittedDensity, col = line.col, lwd = 2, pch = 19)
  points(x = detectedPeaks, y = fittedPeaks, 
         col = peak.col, pch = 19, cex = 1.25)
  points(detectedMins, fittedMins, 
         col = min.col, pch = 19, cex = 1.25)
  text(detectedPeaks, fittedPeaks, col = label.col,
       labels = round(detectedPeaks, 2), pos = 3)
  
}

peakSummaryTable <- function(probe.names.list = NULL, row.id = NULL, 
                             probe.id = "", show.boundary = T,
                             range.start = NULL, range.end = NULL,
                             decimals = 2, var.decimals = 5) {
  if (is.null(range.start)) stop(error = "please specify which rows of the beta matrix the first and last rows of peakSummary correspond to using range.start and range.end.")
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
    if (is.null(probe.names.list)) {
      # Assuming that betas exists in the environment
      probe.id = rownames(betas)[range.start:range.end][row.id]
    } else {
      # Not ideal to pass in a huge character vector
      probe.id = probe.names.list[row.id]
    }
    
  }
  
  title <- paste(title, "for probe", probe.id)
  
  # Row id optional. More useful for debugging
  if (!noRow) {
    title <- paste0(title, " (peakSummary ", row.id, ")")
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

# Updated with new space filtering
# USES BETA MATRIX ONLY: Function to visualize detected peaks compared to histogram
# Assumes that object "betas" exists in the environment
# Assumes that peakSummary does not exist yet
seePeaks <- function(row.id = NULL, probe.id = "", histogram.bins = 50, 
                     density.breaks = 500, label.col = "red", 
                     histogram.col = "black") {
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
    title <- paste0(title, " (Beta matrix row ", row.id, ")")
  } else {
    row.id <- which(rownames(betas) == probe.id)
  }
  
  # fittedDensity <- density(betas[row.id,], from = 0, to = 1, n = density.breaks, 
  #                          adjust = densityAdjust, bw = bandwidth, 
  #                          kernel = kernelType)
  
  calculate <- methylModes(row.data = betas[row.id,])

  print(calculate$detected$propSample)
  
  detectedPeaks <- calculate$detected
  fittedDensity <- calculate$probeDensityEst
  
  #########################
  
  labelHeight <- 0.05*max(fittedDensity$y)
  
  hist(betas[row.id,], breaks = histogram.bins, xlim = c(0,1), 
       probability = TRUE, main = title, col = histogram.col)
  lines(fittedDensity, col = "orchid", lwd = 2, pch = 19)
  points(fittedDensity$x[detectedPeaks$maximaIdx], fittedDensity$y[detectedPeaks$maximaIdx], col = "red", pch = 19, cex = 1.25)
  points(fittedDensity$x[detectedPeaks$leftMinIdx], fittedDensity$y[detectedPeaks$leftMinIdx], col = "blue", pch = 19, cex = 1.25)
  points(fittedDensity$x[detectedPeaks$rightMinIdx], fittedDensity$y[detectedPeaks$rightMinIdx], col = "blue", pch = 19, cex = 1.25)
  text(fittedDensity$x[detectedPeaks$maximaIdx], 
       rep(labelHeight, nrow(detectedPeaks)), col = label.col,
       labels = round(fittedDensity$x[detectedPeaks$maximaIdx], 2), pos = 3)
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
