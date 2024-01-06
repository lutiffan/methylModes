pheno <- read.csv("/nfs/turbo/bakulski1/People/data/ffcw/methyl_data/Feb2022_datafreeze/pd_analytic_freeze1.csv")
colnames(pheno)

# # The whole beta matrix
# betaFileName <- "/nfs/turbo/bakulski1/People/data/ffcw/methyl_data/Feb2022_datafreeze/beta_analytic_freeze1.rds" # file.choose()
# betas <- readRDS(betaFileName) # 5.9 GB

# beta matrix without cross-reactive probes
cleanBetas <- readRDS("cleanedBetas.RDS")
head(colnames(cleanBetas))

# Ensure that the phenotype subject ID ordering matches the beta matrix
sum(colnames(cleanBetas) == pheno$MethID)

#### Filter by child/teen status only ####
isChild <- pheno$childteen == "C"
mean(isChild)
betas <- cleanBetas[, isChild]

isTeen <- pheno$childteen == "T"
mean(isTeen)
betas <- cleanBetas[, isTeen]

# saveRDS(betas, "cleanedBetasChildOnly.RDS")
# teens <- cleanBetas[, !isChild]

#### Filter by sex only ####
isFemale <- pheno$sex == "f"
betas <- cleanBetas[, isFemale]

#### Scratch space ####
betas <- cleanBetas[, isChild & isFemale]
betas <- betas[, isChild & isFemale & pheno$cm1ethrace == 1]
ncol(betas)

#### Filter by chromosome ####
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
relevantLabels <- Locations[rownames(betas),]
yChr <- relevantLabels$chr == "chrY"; sum(yChr)
yChrIdx <- which(yChr)
par(mfrow=c(2,3))
# Take a look at Y chromosome methylation distributions
for (i in 1:6) {
  hist(betas[yChrIdx[i], ], breaks = 50, xlim = c(0,1),
       main = rownames(betas)[yChrIdx][i])
}
dev.off()

# Check out one probe on the X chromosome
xChr <- relevantLabels$chr == "chrX"; sum(xChr)
hist(betas[which(xChr)[1],], breaks = 50)

isAutosomal <- !(yChr | xChr)
# Remove sex chromosomes
cleanBetas <- cleanBetas[isAutosomal,]
# Filter out teens
cleanBetas <- cleanBetas[, isChild]
# saveRDS(cleanBetas, "cleanedBetasChildOnlyAutosomal.RDS")

# Filter out children
cleanBetas <- cleanBetas[, isTeen]
# saveRDS(cleanBetas, "cleanedBetasTeenOnlyAutosomal.RDS")
