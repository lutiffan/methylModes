setwd("/nfs/turbo/bakulski1/People/data/ffcw/methyl_data/Feb2022_datafreeze")
list.files(getwd())
betaFileName <- "beta_analytic_freeze1.rds" # file.choose()
betas <- readRDS(betaFileName) # 5.9 GB

# Check out beta matrix
nrow(betas) # 437588 probes
ncol(betas) # 1684 people
betas[1:3,1:3]

# Mean methylation per probe
probeMeanMeth <- rowMeans(betas)
hist(probeMeanMeth, xlim = c(0,1))

# Mean methylation per person
personMeanMeth <- colMeans(betas)
hist(personMeanMeth, xlim = c(0,1))

# Ignoring phenotypes for now
# phenoBridgeFileName <- "pd_analytic_freeze1.csv" # file.choose()
# phenoBridge <- read.csv(bridgeFileName)
# arrayInfoFileName <- "Methylation_450K_array_batch_information.csv"
# arrayBatchInfo <- read.csv("Methylation_450K_array_batch_information.csv")

# # Try to save space by pre-defining an empty matrix
# nProbes <- nrow(betas)
# nPeople <- ncol(betas)
# # Save an extra column to hold chromosome label
# betaHolder <- matrix(nrow = nProbes, ncol = nPeople + 1)
# # Now fill it (partially) in with the beta matrix
# betaHolder[1:nProbes, 1:nPeople] <- readRDS(betaFileName)

# # Make one histogram
# set.seed(0)
# numProbes <- 9
# samp <- sample(1:nrow(betas), numProbes) # sample some probes
# par(mfrow = c(1,1))
# hist(betas[samp[1],], xlim = c(0,1), breaks = seq(0,1,0.05), probability = T)
# lines(density(betas[samp[1],], kernel = "gaussian"), col = "red")

# Make a lot of histograms - each represents a probe across all people
set.seed(0)
seeds <- sample(1:1000, 50) # Random starting points
# If you get an error about large margins, make the plot window larger
numProbes <- 9
for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  print(paste0("Sampling probes using random seed ",
                                   seeds[i]))
  batch <- sample(1:nrow(betas), numProbes)
  batchProbes <- rownames(betas)[batch]
  par(mfrow = c(3,3)) # Display graphs in a 3x3 grid
  for (j in 1:numProbes) {
    title <- paste0("Beta dist. of probe ", batchProbes[j])
    hist(betas[batch[j],], xlim = c(0,1), breaks = 200, main = title, 
         probability = TRUE)
    lines(density(betas[batch[j],], from = 0, to = 1, adjust = 1), col = "red", lwd = 3)
    print(paste0("Range of data in histogram ", j))
    print(range(betas[batch[j],]))
  }
  browser()
}
dev.off()

######### This code for troubleshooting bad smoothing lines #########
proportionAtEnd <- function(rowNum, left = TRUE, cutoffDist = 0.1, denominator = 1684) {
  if (left) print(sum(betas[rowNum,] < cutoffDist)/denominator)
  else print(sum(betas[rowNum,] > 1 - cutoffDist)/denominator)
}

# Interesting trimodal: row 253061
badIndex <- 73477
hist(betas[badIndex,], xlim = c(0,1), breaks = seq(0,1,0.05), main = title, 
     probability = TRUE)
lines(density(betas[badIndex,], adjust = 1.5, from = 0, to = 1,), col = "red")
#########

# Testing data
testdat <- density(betas[batch[j],], adjust = 1.5, from = 0, to = 1,)
x <- testdat$x
y <- testdat$y

# Try to further smooth the fitted density and count local maxima
# stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset
localMaxima <- function(x, y, half.window = 0.05) {
  require(zoo) # rollapply()
  n <- length(y)
  smoothed <- loess(y ~ x, span = 0.05)$fitted; plot(1:512, smoothed)
  rollingMax <- rollapply(data = zoo(smoothed), 2*half.window + 1, FUN = max, align = "center")
  
  }



if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#View(Islands.UCSC)

# We have the probe names
head(rownames(betas))

# Extract the annotations corresponding to probes in the beta matrix
relevantLabels <- Locations[rownames(betas),]

# Don't do this! It will take forever coercing the beta matrix to characters
# betas[,ncol(betas)] <- relevantLabels$chr

# An ordered vector listing chromosomes: Y, X, and 1-22
chrList <- unique(relevantLabels$chr)
numProbesChr <- numeric(length(chrList))
probeInChrList <- vector(mode = "list", length = length(chrList))

for (i in seq_along(chrList)) {
  print(chrList[i])
  probeInChrList[[i]] <- relevantLabels$chr == chrList[i]
  numProbes <- sum(probeInChrList[[i]])
  numProbesChr[i] <- numProbes
  print(numProbes)
}

barplot(main = "Number of CpGs Measured per Chromosome", height = numProbesChr, 
        names.arg = chrList)

for (i in seq_along(chrList)) {
  hist(rowMeans(betas[probeInChrList[[i]],]))
  browser()
}
