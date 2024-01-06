##### Standard thresholding parameters #####
# Minimum percent of the sample a peak represents
proportionSample = 0.01
# Minimum space between peaks on the x axis
personalSpace = 0.10
# Between two adjacent peaks, this is the maximum ratio of the minima height to 
# the shorter of two peaks with nearly equal or equal height
maxDipRatio = 0.75
minDipRatio = 0.50

##### Advanced parameters for smoothing #####
bwOptions <- c(NA, "sheatherJones")
kernelType = "gaussian" # "This must partially match one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine", with default "gaussian", and may be abbreviated to a unique prefix (single letter)."
bandwidthType = bwOptions[1] # options: NA, "sheatherJones"
numBreaks = 500 # The "n" parameter for density(): "the number of equally spaced points at which the density is to be estimated..." You need at least 500 for a decent smooth curve
densityAdjust = 1
pushToZero = 1e-6 # Get rid of minuscule floating point numbers  created by density() that really should be zero

# Applicable when using subset of beta matrix
# rangeStart = 1
# rangeEnd = 401545
# totalRows = rangeEnd - rangeStart + 1