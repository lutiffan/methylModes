# Function to return indices of local min/maxima

# # Can use simulated data for debugging
# fitted <- probeDensityEst$y
# zeroThreshold = 0.01

# Assume input data is not all slope == 0
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
  
  data.frame("maximaIdx" = slopeZero[whereMax], "leftMinIdx" = c(1, slopeZero[leftMin]), "rightMinIdx" = slopeZero[rightMin])
}
