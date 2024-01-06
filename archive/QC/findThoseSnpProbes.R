# Needs 75+ GiB memory to load in meth object
# Need 130+ GB memory to run QC steps
library(dplyr)
library(minfi)
source("visualizeBetaPeaks.R")

# snpbeta <- readRDS("~/Peak_detection_old/snpbeta.RDS")
# 
# annotations <- read.csv(file = "/nfs/turbo/bakulski1/People/data/ffcw/methyl_data/Feb2022_datafreeze/pd_analytic_freeze1.csv") 
# cleansnpbeta <- snpbeta[,colnames(snpbeta) %in% annotations$MethID]
# 
# numNAs <- numeric(nrow(cleansnpbeta))
# for ( i in 1:nrow(cleansnpbeta)) {
#   numNAs[i] <- sum(is.na(cleansnpbeta[i,]))
# }
# 
# # Why do we still have NAs? Cut them for now
# cleansnpbeta <- cleansnpbeta[numNAs == 0,]
# 
# # Filter out one of the age groups to avoid near-duplicate data
# pheno <- read.csv("/nfs/turbo/bakulski1/People/data/ffcw/methyl_data/Feb2022_datafreeze/pd_analytic_freeze1.csv")
# isChild <- pheno$childteen == "C"
# cleansnpbeta <- cleansnpbeta[, isChild]

# saveRDS(cleansnpbeta, "cleansnpbeta.RDS")
cleansnpbeta <- readRDS("cleansnpbeta.RDS")

ghSNPs <- gaphunter(object = cleansnpbeta, verbose = T)

ghCounts <- table(ghSNPs$proberesults$Groups)
barplot(ghCounts, main = "Gaphunter Group Counts Among 44 SNP Probes")

betas <- cleansnpbeta
for (i in 1:nrow(cleansnpbeta)) {
  seePeaks(row.id = i)
  legend("top", as.character(ghSNPs$proberesults$Groups[i]), pch=1, title= "Gaphunter group count: ", inset = .05)
  browser()
}
