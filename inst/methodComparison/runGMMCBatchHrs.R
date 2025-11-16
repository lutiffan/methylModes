library(mclust)
library(data.table)
library(qs)
library(parallelly)
library(doParallel)
library(iterators)

# Test command
# sbatch --export=begin=1,end=10000 /nfs/turbo/bakulski1/People/lutiffan/hrs_analyses/gaussianMixtureModelClustering/runGMMCBatchHrs.sh

#### Setup ####
beginIdx = Sys.getenv("begin")
endIdx = Sys.getenv("end")

# Load beta matrix
betas <- qread("/nfs/turbo/bakulski1/Datasets/HRS/DNA_methylation/Cleaned_Apr2024/HRS_DNAm.qs")
betas <- betas[beginIdx:endIdx,]

startTime <- Sys.time()

# Make new cluster
cl <- parallelly::makeClusterPSOCK(availableCores(omit = 1), autoStop = TRUE)

# Register parallel backend
doParallel::registerDoParallel(cl)

template <- data.table("probeId" = character(1),
                       "optimalG" = numeric(1),
                       "mixingProportions" = vector(mode = "list", 
                                                length = 1),
                       "sampleProportions" = vector(mode = "list", length = 1),
                       "means" = vector(mode = "list", length = 1),
                       "sigmasq" = vector(mode = "list", length = 1))

gmmcResults <- foreach::foreach(probe = iter(betas, by = "row"), 
                                .combine = "rbind",
                                .packages = c("foreach", "data.table", "mclust"),
                                .export = c("betas")) %dopar% {
  fit <- Mclust(data = probe, G = 1:3, verbose = F)
  # z is a matrix whose [i,k]th entry is the prob. that obs. i belongs to class k
  # The dimensions of fit$z are n x (estimated # clusters)
  mixingProportions <- colSums(fit$z)/sum(fit$z)
  # length(fit$classification) = n
  sampleProportions <- as.numeric(table(fit$classification))/fit$n
  
  template[1,] <- list("probeId" = rownames(probe),
       "optimalG" = fit$G,
       "mixingProportions" = mixingProportions,
       "sampleProportions" = sampleProportions,
       "means" = fit$parameters$mean[!is.na(fit$parameters$mean)],
       "sigmasq" = fit$parameters$variance$sigmasq[!is.na(fit$parameters$variance$sigmasq)])
  template
  }
stopCluster(cl)
elapsed = Sys.time() - startTime

# # Create container for results
# gmmcResults <- data.table("probeId" = character(totalRows),
#                           "optimalG" = numeric(totalRows), 
#                           "mixingProportions" = vector(mode = "list", 
#                                                 length = totalRows),
#                           "sampleProportions" = vector(mode = "list",
#                                                 length = totalRows),
#                           "means" = vector(mode = "list", 
#                                            length = totalRows),
#                           "sigmasq" = vector(mode = "list", 
#                                              length = totalRows),
#                           "classification" = vector(mode = "list", 
#                                                     length = totalRows))
# 
# #### Main for-loop ####
# 
# for (i in 1:totalRows) {
#   
  # fit <- Mclust(data = betas[i,], G = 1:3, verbose = F)
  # 
  # # z is a matrix whose [i,k]th entry is the prob. that obs. i belongs to class k
  # # The dimensions of fit$z are n x (estimated # clusters)
  # mixingProportions <- colSums(fit$z)/sum(fit$z)
  # # length(fit$classification) = n
  # sampleProportions <- as.numeric(table(fit$classification))/fit$n
  # 
  # gmmcResults[i,] <- list("probeId" = rownames(betas)[i],
  #                         "optimalG" = fit$G,
  #                         "mixingProportions" = mixingProportions,
  #                         "sampleProportions" = sampleProportions,
  #                         "means" = fit$parameters$mean[!is.na(fit$parameters$mean)],
  #                         "sigmasq" = fit$parameters$variance$sigmasq[!is.na(fit$parameters$variance$sigmasq)],
  #                         "classification" = fit$classification)
# }

#### Logging and saving ####
output <- paste0("/nfs/turbo/bakulski1/People/lutiffan/hrs_analyses/gaussianMixtureModelClustering/results/", 
                 "begin_", format(beginIdx, scientific = F),
                 "_end_", format(endIdx, scientific = F),
                 "_date_", Sys.Date(), ".RDS")
# Keep track of runtimes
sink("/nfs/turbo/bakulski1/People/lutiffan/hrs_analyses/gaussianMixtureModelClustering/results/runtimesGMMC.txt", append = TRUE)

print(output)
elapsed

sink()

saveRDS(gmmcResults, file = output)
