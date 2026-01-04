#' Detect Local Minima and Maxima in a Density Function
#' 
#' This function identifies the indices of local minima and maxima in a density function,
#' typically used for peak detection in DNA methylation data. It handles edge cases and
#' ensures proper boundary conditions for the density estimation.
#' 
#' @param methylModesResults A data.table returned by methylModesBatch()
#' 
#' @return A modified version of the input data table with these changes:
#' \enumerate{
#'   \item Within each row, columns peakLocations, proportionSample, and 
#'   peakVariance (which are objects of class "list") are sorted in descending 
#'   order by proportionSample
#'   \item Columns are appended to the end of the data frame. For peakLocations, 
#'   proportionSample, and peakVariance, the values listed within these columns
#'   are split across new columns containing one value for each peak. The 
#'   number of new columns is equal to 3 * the maximum number of peaks found for
#'   any probe. These columns are named with the following convention:
#'   peakLocation_i, proportionSample_i, peakVariance_i, where i ranges from 1 
#'   to the maximum number of peaks
#' }
#' 
#' @examples
#' # Obtain MethylModes results from the toy dataset
#' # Note that this will take a minute or so
#' test <- methylModesBatch(methylModes::toyMultimodalData)
#' testExpanded <- expandPeakDataCols(test)
#' # Sort by second-largest peak as a measure of confidence in multimodality
#' testExpanded <- data.table::setorder(testExpanded, -proportionSample_2, na.last = TRUE)
#' 
#' @export

expandPeakDataCols <- function(methylModesResults) {
  # Ensure data.table is used
  setDT(methylModesResults)

  # Count number of rows where methylModesResults$numPeaks is NA
  numRowsWithNA <- sum(is.na(methylModesResults$numPeaks))
  # Print a message to the console stating the number of rows removed
  if (numRowsWithNA > 0) {
    message(paste0("Removed ", numRowsWithNA, " rows where numPeaks is NA"))
  }
  
  # Filter out rows where methylModesResults$numPeaks is NA
  methylModesResults <- methylModesResults[!is.na(methylModesResults$numPeaks),]
  # Print a message to the console stating the number of rows removed
  
  # Validate peak number
  maxPeaks <- max(methylModesResults$numPeaks)
  
  # Sort peaks by proportionSample (largest to smallest) for each probe
  methylModesResults[, c("peakLocations", "proportionSample", "peakVariance") := {
    # Get order indices sorted by proportionSample descending
    order_idx <- lapply(proportionSample, function(x) order(-x))
    
    # Sort each list according to the order
    list(
      mapply(function(loc, idx) loc[idx], peakLocations, order_idx, SIMPLIFY = FALSE),
      mapply(function(prop, idx) prop[idx], proportionSample, order_idx, SIMPLIFY = FALSE),
      mapply(function(var, idx) var[idx], peakVariance, order_idx, SIMPLIFY = FALSE)
    )
  }]
  
  # Expand list columns into long format first
  expanded <- methylModesResults[, .(
    peak = seq_len(numPeaks),
    peakLocation = unlist(peakLocations),
    proportionSample = unlist(proportionSample),
    peakVariance = unlist(peakVariance)
  ), by = probeName]
  
  # Reshape from long to wide format
  wide <- dcast(
    expanded,
    probeName ~ peak,
    value.var = c("peakLocation", "proportionSample", "peakVariance")
  )
  
  # Merge back with original
  result <- merge(methylModesResults, wide, by.x = "probeName", by.y = "probeName", all.x = TRUE)
  
  return(result[])
}

# # Tooooo slow but the logic is easy to follow for a reader
# expandPeakDataCols <- function(methylModesResults, maxPeakNumber = NA) {
#   if ( (maxPeakNumber < 1) | (maxPeakNumber > max(methylModesResults$numPeaks)) ) {
#     stop("Peak number must be at least one and at most maximum number of peaks identified in this set of results.")
#   }
#   if(is.na(maxPeakNumber)) maxPeakNumber = max(methylModesResults$numPeaks)

#   for (c in 1:maxPeakNumber) {
#     methylModesResults[, paste0("peakLocations", c) := as.numeric(NA)]
#     methylModesResults[, paste0("proportionSample", c) := as.numeric(NA)]
#     methylModesResults[, paste0("peakVariance", c) := as.numeric(NA)]
#   }

#   for (i in 1:nrow(methylModesResults)) {
#     locations <- unlist(methylModesResults$peakLocations[[i]])
#     proportions <- unlist(methylModesResults$proportionSample[[i]])
#     variances <- unlist(methylModesResults$peakVariance[[i]])

#     for (p in 1:methylModesResults$numPeaks[i]) {
#       methylModesResults[i, (paste0("peakLocations", p)) := locations[p]]
#       methylModesResults[i, (paste0("proportionSample", p)) := proportions[p]]
#       methylModesResults[i, (paste0("peakVariance", p)) := variances[p]]
#     }
#   }
# }