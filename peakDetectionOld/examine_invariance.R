##### Data reading and writing #####

setwd("/nfs/turbo/bakulski1/People/data/ffcw/methyl_data/Feb2022_datafreeze")
list.files(getwd())
betaFileName <- "beta_analytic_freeze1.rds" # file.choose()
betas <- readRDS(betaFileName) # 5.9 GB
library(data.table)

betas <- as.data.table(betas, keep.rownames = TRUE) # A little smaller

# tbetas <- transpose(betas) # slightly better

# betas$id <- 1:nrow(betas)
# tbeta <- dcast(melt(betas, id.vars = "id"), variable ~ id) # too slow



##### Mean methylation per probe #####

probeMeanMeth <- rowMeans(betas)
hist(probeMeanMeth, xlim = c(0,1), breaks = 300, probability = TRUE)
lines(density(probeMeanMeth, from = 0, to = 1), col = "red")
isHypo <- probeMeanMeth < 0.3
isHyper <- probeMeanMeth > 0.7
mean(isHypo) # 0.4240816
mean(isHyper) # 0.4734499 # But the variance is larger at this end
var(probeMeanMeth[isHypo]) # 0.004459489
var(probeMeanMeth[isHyper]) # 0.005419683

# For reference and curiosity
format(object.size(probeMeanMeth), units = "MB")

threshold <- 0.1

# Proportion of largely unmethylated probes
# All betas across all subjects
sum(betas < threshold)/length(betas) # for threshold = 0.1: 0.3708118
# Compare to averaging probes across all subjects
sum(probeMeanMeth < threshold)/length(probeMeanMeth) # for threshold = 0.1: 0.3625922

# Proportion of largely methylated probes
# All betas across all subjects
sum(betas > 1 - threshold)/length(betas) # for threshold = 0.1 -> 0.3277098
# Compare to averaging probes across all subjects
sum(probeMeanMeth > 1 - threshold)/length(probeMeanMeth) # for threshold = 0.1: 0.3172733

# Check out some interesting probes
sum(betas[211484,] < 0.1)/ncol(betas) # Trimodal, tall peaks
sum(betas[179171,] < 0.1)/ncol(betas) # Trimodal, short peaks
sum(betas[rownames(betas) == 'cg07419801',] < 0.1)/ncol(betas) 



##### Proportion of probes with close to invariant methylation #####

# methylCategory <- function(probe, denominator = 1684) {
#   categories <- numeric(2)
#   categories[1] <- mean(probe < 0.1)
#   categories[2] <- mean(probe > 0.9)
#   
#   categories
# }

betaSummary <- matrix(nrow = nrow(betas), ncol = 3)
colnames(betaSummary) <- c("hypo", "mid", "hyper")

threshold <- 0.3

for (i in 1:nrow(betas)) {
  hypo <- mean(betas[i,] < threshold)
  hyper <- mean(betas[i,] > 1 - threshold)
  mid <- 1 - hypo - hyper
  betaSummary[i,] <- c(hypo, mid, hyper)
}

# Examine hypo and hyper-methylated sites
# for threshold = 0.1, 0.2653501 are cleanly binarized to hypo and hyper
# for threshold = 0.3, 0.6027268 are in the hypo- and hyper- methyl categories
mean(betaSummary[,2] == 0) 

# Proportion of probes with samples in all three categories
# for threshold = 0.3, 0.1283765
mean(betaSummary[,1] > 0 & betaSummary[,2] > 0 & betaSummary[,3] > 0)

# Proportion of probes that are only hypomethylated
# for threshold = 0.3, 0.3450643
hypoOnly <- betaSummary[,1] == 1
mean(hypoOnly)

# Proportion of probes that are just in the middle
# Very rare. for threshold = 0.3, 0.001092352
mean(betaSummary[,2] == 1) 

# Proportion of probes that are only hypermethylated
# for threshold = 0.3, 0.2523584
hyperOnly <- betaSummary[,3] == 1
mean(hyperOnly)

# Experiment with data.table functionality
betaSummary2 <- betas[lapply(.SD), ]
# isInvariant <- ifelse( (betas > (1 - threshold)) | (betas < threshold), 
#                        TRUE, FALSE)
