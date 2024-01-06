library(mclust)
library(data.table)

#### Command line parameters ####
# Specify range of row numbers of beta matrix
dataset <- Sys.getenv("data")
rangeStart = as.numeric(Sys.getenv("start"))
rangeEnd = as.numeric(Sys.getenv("end"))
selection = Sys.getenv("selection") # Selection criteria

# # Testing
# dataset = "child"
# rangeStart = 1
# rangeEnd = 10
# selection = "BIC"

selectionMethods <- c("BIC", "LR")
if (!(selection %in% selectionMethods)) {
  stop(simpleError(paste("Choose one of the model selection methods:", selectionMethods)))
}

#### Setup ####
# Load beta matrix
if (dataset == "child") {
  betas <- readRDS("/home/lutiffan/betaMatrix/cleanedBetasChildOnlyAutosomal.RDS")
} else if (dataset == "teen") {
  betas <- readRDS("/home/lutiffan/betaMatrix/cleanedBetasTeenOnlyAutosomal.RDS")
} else {
  stop(simpleError("invalid dataset name."))
}

# Sample size
n <- ncol(betas)

# Subset it to the range we want to cover in a given batch
betas <- betas[rangeStart:rangeEnd,]

# Create container for results
totalRows <- nrow(betas)
gmmcResults <- data.table("probeId" = character(totalRows),
                          "optimalG" = numeric(totalRows), 
                          "mixingProportions" = vector(mode = "list", 
                                                length = totalRows),
                          "sampleProportions" = vector(mode = "list",
                                                length = totalRows),
                          "means" = vector(mode = "list", 
                                           length = totalRows),
                          "sigmasq" = vector(mode = "list", 
                                             length = totalRows),
                          "classification" = vector(mode = "list", 
                                                    length = totalRows))

#### Main for-loop ####
startTime <- Sys.time()

if (selection == "BIC") {
  for (i in 1:totalRows) {
    
    fit <- Mclust(data = betas[i,], G = 1:3, verbose = F)
    
    mixingProportions <- colSums(fit$z)/sum(fit$z)
    sampleProportions <- as.numeric(table(fit$classification))/n
    
    gmmcResults[i,] <- list("probeId" = rownames(betas)[i],
                            "optimalG" = fit$G,
                            "mixingProportions" = mixingProportions,
                            "sampleProportions" = sampleProportions,
                            "means" = fit$parameters$mean[!is.na(fit$parameters$mean)],
                            "sigmasq" = fit$parameters$variance$sigmasq[!is.na(fit$parameters$variance$sigmasq)],
                            "classification" = fit$classification)
    
    # # Sometimes the clusters are ridiculously small.
    # keepCluster <- as.numeric(which(table(fit$classification)/ncol(betas) > proportionSample))
    # 
    # gmmcResults[i,] <- list("probeId" = rownames(betas)[i],
    #                         "optimalG" = length(keepCluster),
    #                         "groupNames" = keepCluster,
    #                         "means" = fit$parameters$mean[keepCluster],
    #                         "sigmasq" = fit$parameters$variance$sigmasq[keepCluster],
    #                         "classification" = fit$classification,
    #                         "G" = fit$G)
  }
} else {
  stop(simpleError("Likelihood ratio testing under construction."))
}

#### Logging and saving ####
output <- paste0("/home/lutiffan/GMMCResults/", dataset, "_", selection, "_", 
                 format(rangeStart, scientific = F),"_", 
                 format(rangeEnd, scientific = F), ".RDS")
# Keep track of runtimes
sink("/home/lutiffan/GMMCResults/runtimesGMMC.txt", append = TRUE)
Sys.time()
print(output)
Sys.time() - startTime
sink()

saveRDS(gmmcResults, file = output)