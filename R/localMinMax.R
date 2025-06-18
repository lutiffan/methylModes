#' Detect Local Minima and Maxima in a Density Function
#' 
#' This function identifies the indices of local minima and maxima in a density function,
#' typically used for peak detection in DNA methylation data. It handles edge cases and
#' ensures proper boundary conditions for the density estimation.
#' 
#' @param fitted A numeric vector representing the density function values.
#' @param zeroThreshold A numeric value (default: 1e-6) below which values are considered zero.
#' 
#' @return A data frame with three columns:
#' \describe{
#'   \item{maximaIdx}{Indices of local maxima in the input vector}
#'   \item{leftMinIdx}{Indices of left minima for each maximum}
#'   \item{rightMinIdx}{Indices of right minima for each maximum}
#' }
#' 
#' @details
#' The function processes the input vector in several steps:
#' \enumerate{
#'   \item Values below zeroThreshold are set to zero
#'   \item Slope changes are detected using diff()
#'   \item Local maxima and minima are identified based on slope changes
#'   \item Edge cases are handled (e.g., flat regions at start/end)
#' }
#' 
#' @examples
#' # Create a simple density function with two peaks
#' x <- seq(0, 1, length.out = 100)
#' y <- dnorm(x, mean = 0.3, sd = 0.1) + dnorm(x, mean = 0.7, sd = 0.1)
#' 
#' # Find local minima and maxima
#' peaks <- localMinMax(y)
#' 
#' # Plot the results
#' plot(x, y, type = "l", main = "Density Function with Detected Peaks")
#' points(x[peaks$maximaIdx], y[peaks$maximaIdx], col = "red", pch = 16, cex = 1.5)
#' points(x[peaks$leftMinIdx], y[peaks$leftMinIdx], col = "blue", pch = 16)
#' points(x[peaks$rightMinIdx], y[peaks$rightMinIdx], col = "blue", pch = 16)
#' legend("topright", legend = c("Maxima", "Minima"), 
#'        col = c("red", "blue"), pch = 16)
#' 
#' # Print the detected peaks
#' cat("Detected", nrow(peaks), "peaks:\n")
#' for (i in 1:nrow(peaks)) {
#'   cat("Peak", i, ": Maximum at index", peaks$maximaIdx[i], 
#'       "(x =", round(x[peaks$maximaIdx[i]], 3), ")\n")
#' }
#' 

localMinMax <- function(fitted, zeroThreshold = 1e-6) {
  
  fitted <- ifelse(fitted < zeroThreshold, 0, fitted)
  
  # Cases:
  # 1. If line increasing at start, plusMinus will be TRUE TRUE...
  # until we reach the first stationary point to the right
  
  # 2. If line decreasing/flat at start, plusMinus will be TRUE FALSE ... 
  # detect stationary point at index == 1
  # When slope is no longer zero, we will incorrectly detect a stationary point
  
  # this code cannot distinguish between negative slope and slope of zero
  plusMinus <- diff(c(-Inf, fitted)) > 0L
  
  # Find indices where the sign of the slope changes
  slopeZero <- cumsum(rle(plusMinus)$lengths); slopeZero
  
  if (length(slopeZero) == 1) {
    return(list("maximaIdx" = slopeZero, "minimaIdx" = c(1, 500)))
  }
  
  # 1. First slope sign change is a max
  # 2. First (artificial) slope sign change is a max and is located at index == 1
  # 3. First (artificial) slope sign change is just artificial. Second sign change is an inflection point
  whereMax <- seq.int(1L, length(slopeZero), 2L)
  whereMin <- (1:length(slopeZero))[-whereMax] 
  
  # Case 3: Line starts out flat
  if (fitted[[1]] == fitted[[2]]) { 
    # Skip artificial first slope sign change (since we start from -Inf)
    whereMax <- whereMax[-1] 
    
    # Skip inflection point
    whereMin <- whereMin[-1]
  }
  
  leftMin <- whereMax - 1
  # We will force first minima to be at zero so detected max/min spans the whole 0 to 1 range
  leftMin <- leftMin[-1]
  
  # # Handle edge case where first max is at first index
  # if (leftMin[1] == 0) leftMin[1] = 1
  
  
  rightMin <- whereMax + 1
  # Handle edge case where last max is at last index
  if (rightMin[length(rightMin)] > length(slopeZero)) {
    rightMin[length(rightMin)] = length(slopeZero)
  } 
  
  data.frame("maximaIdx" = slopeZero[whereMax], 
             "leftMinIdx" = c(1, slopeZero[leftMin]), 
             "rightMinIdx" = slopeZero[rightMin])
}
