library(data.table)
library(dplyr)

manifestEpic <- fread("/nfs/turbo/bakulski1/People/lutiffan/dataCleaningIlluminaManifests/MethylationEPIC_v-1-0_B4.csv", skip = 7)
neededColumnsEpic <- c("IlmnID", "CHR", "Relation_to_UCSC_CpG_Island", "MAPINFO", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")
epicInfo1 <- manifestEpic[1:865918, ..neededColumnsEpic]

# Convert character strings to numeric lists
snpDistanceNumericAll <- vector(mode = "list", length = nrow(epicInfo1))
# mafNumericAll <- vector(mode = "list", length = nrow(epicInfo1))
# numSnps <- numeric(nrow(epicInfo1))

for (i in seq_along(snpDistanceNumericAll)) {
  if (epicInfo1$SNP_DISTANCE[i] == "") {
    snpDistanceNumericAll[i] <- NA
    # mafNumericAll[i] <- NA
    next
  }
  snpDists <- base::strsplit(epicInfo1$SNP_DISTANCE[i], split = "[;]")
  snpDistanceNumericAll[[i]] <- as.numeric(unlist(snpDists))
  
  # mafs <- base::strsplit(epicInfo1$SNP_MinorAlleleFrequency[i], split = "[;]")
  # mafNumericAll[[i]] <- as.numeric(unlist(mafs))
  
  # numSnps[i] <- length(snpDistanceNumericAll[[i]])
}

epicInfo1$snpDistanceLists <- snpDistanceNumericAll

# Create three columns counting number of SNPs falling into the following categories:
# 1. 0-1
# 2. 2-10
# 3. 11-50
distanceCategoryCount <- function(x, min, max) {
    x <- unlist(x)
    return(sum(x >= min & x <= max))
}
epicInfo1$snpDistance0_1 <- sapply(snpDistanceNumericAll, distanceCategoryCount, min = 0, max = 1)
epicInfo1$snpDistance2_10 <- sapply(snpDistanceNumericAll, distanceCategoryCount, min = 2, max = 10)
epicInfo1$snpDistance11_50 <- sapply(snpDistanceNumericAll, distanceCategoryCount, min = 11, max = 50)

IlluminaManifestEPIC <- epicInfo1 %>%
  select(-SNP_DISTANCE) %>%
  rename(SNP_DISTANCE = snpDistanceLists) 

save(IlluminaManifestEPIC, file = "IlluminaManifestEPIC.rda")