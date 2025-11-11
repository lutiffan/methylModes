library(data.table)
library(dplyr)
source("closestIsland.R")

manifestEpic2 <- fread("/nfs/turbo/bakulski1/People/lutiffan/dataCleaningIlluminaManifests/EPIC-8v2-0_A2.csv", skip = 7)
neededColumnsEpic2 <- c("Name", "CHR", "Relation_to_UCSC_CpG_Island", "UCSC_CpG_Islands_Name", "MAPINFO", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")
# sum(manifestEpic2$Probe_Type == "cg") - sum(manifestEpic2$IlmnID %like% "cg" & manifestEpic2$CHR > 0) # 6881
# sum(manifestEpic2$Probe_Type == "cg" & manifestEpic2$CHR ==0) # 6881

epicInfo2 <- manifestEpic2[manifestEpic2$IlmnID %like% "cg" & manifestEpic2$CHR > 0, ..neededColumnsEpic2]
rm(manifestEpic2)

# # Get closest island for every CpG
# closestIslandIndex <- numeric(nrow(epicInfo2))
# for (i in seq_along(closestIslandIndex)) {
#   closestIslandIndex[i] <- closestIsland(epicInfo2$MAPINFO[i], epicInfo2$UCSC_CpG_Islands_Name[i])
# }

# Seems that all of the islands are sorted in order of nearest to farthest from the CpG
# So we can just take the first island for each CpG
epicInfo2$closest_CpG_Island_Feature <- sapply(epicInfo2$Relation_to_UCSC_CpG_Island, function(x) strsplit(x, split = ";")[[1]][1])
epicInfo2$UCSC_CpG_Islands_Name[epicInfo2$UCSC_CpG_Islands_Name == ""] <- "OpenSea"

# Convert character strings to numeric lists
snpDistanceNumericAll <- vector(mode = "list", length = nrow(epicInfo2))
# mafNumericAll <- vector(mode = "list", length = nrow(epicInfo2))
# numSnps <- numeric(nrow(epicInfo2))

for (i in seq_along(snpDistanceNumericAll)) {
  if (epicInfo2$SNP_DISTANCE[i] == "") {
    snpDistanceNumericAll[i] <- NA
    # mafNumericAll[i] <- NA
    next
  }
  snpDists <- base::strsplit(epicInfo2$SNP_DISTANCE[i], split = "[;]")
  snpDistanceNumericAll[[i]] <- as.numeric(unlist(snpDists))
  
  # mafs <- base::strsplit(epicInfo2$SNP_MinorAlleleFrequency[i], split = "[;]")
  # mafNumericAll[[i]] <- as.numeric(unlist(mafs))
  
  # numSnps[i] <- length(snpDistanceNumericAll[[i]])
}

epicInfo2$snpDistanceLists <- snpDistanceNumericAll
# epicInfo2$snpMinorAlleleFreq <- mafNumericAll
# epicInfo2$numSnps <- numSnps

# # Extract associated info for the nearest SNP
# minDist <- numeric(nrow(epicInfo2))
# minDistMAF <- numeric(nrow(epicInfo2))
# for (i in seq_along(minDist)) {
#   if (epicInfo2$SNP_DISTANCE[i] == "") {
#     minDist[i] <- NA
#     minDistMAF[i] <- NA
#     next
#   }
  
#   closest <- which.min(epicInfo2$snpDistanceLists[[i]])
#   minDist[i] <- epicInfo2$snpDistanceLists[[i]][closest]
#   minDistMAF[i] <- epicInfo2$snpMinorAlleleFreq[[i]][closest]
# }

# epicInfo2$minDist <- minDist
# epicInfo2$minDistMAF <- minDistMAF

# Create three columns counting number of SNPs falling into the following categories:
# 1. 0-1
# 2. 2-10
# 3. 11-50
distanceCategoryCount <- function(x, min, max) {
    x <- unlist(x)
    return(sum(x >= min & x <= max))
}
epicInfo2$snpDistance0_1 <- sapply(snpDistanceNumericAll, distanceCategoryCount, min = 0, max = 1)
epicInfo2$snpDistance2_10 <- sapply(snpDistanceNumericAll, distanceCategoryCount, min = 2, max = 10)
epicInfo2$snpDistance11_50 <- sapply(snpDistanceNumericAll, distanceCategoryCount, min = 11, max = 50)


IlluminaManifestEPICv2 <- epicInfo2 %>% select(Name, 
CHR, 
closest_CpG_Island_Feature, 
MAPINFO, 
SNP_ID, 
snpDistanceLists,
snpDistance0_1,
snpDistance2_10,
snpDistance11_50) %>%
  rename(IlmnID = Name,
    Relation_to_UCSC_CpG_Island = closest_CpG_Island_Feature,
    SNP_DISTANCE = snpDistanceLists)

save(IlluminaManifestEPICv2, file = "IlluminaManifestEPICv2.rda")
