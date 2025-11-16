postProcessing <- function(gmmcResults, 
                           OUT_CUTT = NA, 
                           MAX_STD = 0.1, 
                           MIN_MEAN_DIFF = 0.2) {
  gmmcResults[, `:=`(
    # "First, the largest cluster cannot be too big"
    # This equals T if largest cluster contains too high of a proportion of the sample
    bigCluster = sapply(sampleProportions, function(x) max(unlist(x)) > (1 - OUT_CUTT)),
    # Second, samples within each cluster should have small variance
    # This equals T if any of the clusters have too large of a variance
    bigVariance = sapply(sigmasq, function(x) any(unlist(x) > MAX_STD)),
    # Finally, we require that cluster centers should be separable from each other
    # This equals T if any two clusters have a distance less than allowed
    tooClose = sapply(means, function(x) any(diff(unlist(x)) < MIN_MEAN_DIFF))
  )]
  gmmcResults[!(gmmcResults$bigCluster | gmmcResults$bigVariance | gmmcResults$tooClose),]
}