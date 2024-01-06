library(data.table)
library(dplyr)
source("/home/lutiffan/peakDetectionScripts/fillPeakSummary.R") # sources peakPicker.R, which sources localMinMax.R
# source("iterateOverProbes.R")

# Load beta matrix
# betaFileName <- "/nfs/turbo/bakulski1/People/data/ffcw/methyl_data/Feb2022_datafreeze/beta_analytic_freeze1.rds" # file.choose()
# betas <- readRDS(betaFileName) # 5.9 GB

# # Set parameters manually for testing
# dataset = "child"
# rangeStart = 3000
# rangeEnd = 3300
# propStart = 0.025
# propEnd = 0.025
# spaceStart = 0.025
# spaceEnd = 0.025
# spacing = 0.025

##### Adjustable parameters #####
# Choose child or teen data
dataset <- Sys.getenv("data")
# Specify range of row numbers of beta matrix
rangeStart = as.numeric(Sys.getenv("start"))
rangeEnd = as.numeric(Sys.getenv("end"))
# Specify a value or a range of values for proportionSample
propStart = as.numeric(Sys.getenv("prstart")) # 0.025
# To run one value of proportionSample, set propEnd to same value as propStart
propEnd = as.numeric(Sys.getenv("prend")) # 0.1
# Specify a value or a range of values for personalSpace
spaceStart = as.numeric(Sys.getenv("spstart")) # 0.025
# To run one value of personalSpace, set propEnd to same value as personalSpace
spaceEnd = as.numeric(Sys.getenv("spend")) # 0.1

maxDipRatio = as.numeric(Sys.getenv("dr"))
# For an evenly spaced sequence of parameter values, specify spacing size
spacing = as.numeric(Sys.getenv("spacing")) # 0.025
destination = Sys.getenv("destination")

if (!(destination %in% c("sensitivityAnalysis", "peakDetectionResults"))) {
  stop(simpleError("Invalid output directory."))
}

if (dataset == "child") {
  betas <- readRDS("/home/lutiffan/betaMatrix/cleanedBetasChildOnlyAutosomal.RDS")
} else if (dataset == "teen") {
  betas <- readRDS("/home/lutiffan/betaMatrix/cleanedBetasTeenOnlyAutosomal.RDS")
}


if (propStart == propEnd) {
  propValues <- propStart
} else {
  propValues <- seq(from = propStart, to = propEnd, by = spacing)
}

if (spaceStart == spaceEnd) {
  spaceValues <- spaceStart
} else {
  spaceValues <- seq(from = spaceStart, to = spaceEnd, by = spacing)
}

# Advanced parameters for smoothing
bandwidthType = NA # options: NA, "SJ"
kernelType = "gaussian"
numBreaks = 500 # The "n" parameter for density(): "the number of equally spaced points at which the density is to be estimated..." You need at least 500 for a decent smooth curve
densityAdjust = 1
pushToZero = 1/ncol(betas) # 1e-10 # Get rid of minuscule floating point numbers created by density() that really should be zero

#################################

totalRows = rangeEnd - rangeStart + 1

for (p in seq_along(propValues)) {
  proportionSample = propValues[p]
  
  # Create empty container for probe-level summaries
  peakSummary <- data.table("numPeaks" = numeric(totalRows), 
                            "meanBeta" = numeric(totalRows),
                            "peakLocations" = vector(mode = "list", 
                                                     length = totalRows),
                            "leftMin" = vector(mode = "list", 
                                               length = totalRows),
                            "rightMin" = vector(mode = "list", 
                                                length = totalRows),
                            "proportionSample" = vector(mode = "list", 
                                                        length = totalRows),
                            "peakVariance" = vector(mode = "list", 
                                                    length = totalRows),
                            "gapFound" = logical(totalRows))
  
  for (s in seq_along(spaceValues)) {
    personalSpace = spaceValues[s]
    
    startTime <- Sys.time()
    peakSummary <- fillPeakSummary()
    endTime <- Sys.time() - startTime
    
    # Output filename
    output <- paste0("/home/lutiffan/", destination, "/peakSummary_", 
                     dataset, "_autosomal_",
                     gsub('\\.', '_', proportionSample), "_space_",
                     gsub('\\.', '_', personalSpace), "_start_",
                     format(rangeStart, scientific = F), "_end_", 
                     format(rangeEnd, scientific = F), ".RDS")
    
    # output <- paste0("sensitivityAnalysis/test_peakSummary_prop_",
    #                  gsub('\\.', '_', proportionSample), "_space_",
    #                  gsub('\\.', '_', personalSpace), "_start_",
    #                  rangeStart, "_end_", rangeEnd, ".RDS")
    
    # Keep track of how long it took to fill peakSummary
    sink(paste0(destination, "/runtimes.txt"), append = TRUE)
    print(output)
    print(endTime)
    sink()
    
    saveRDS(peakSummary, file = output)
  }
}
