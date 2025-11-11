library(data.table)
library(dplyr)
# Line 485586 is the start of the list of controls
manifest450k <- fread("/nfs/turbo/bakulski1/People/lutiffan/dataCleaningIlluminaManifests/humanmethylation450_15017482_v1-2.csv", skip = 7)

# Check against bioconductor 450k annotation library
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
Locations[(rownames(Locations) %in% c("cg00035864", "cg24687794", "cg01751584")),]

View(manifest450k[manifest450k$IlmnID %in% c("cg00035864", "cg24687794", "cg01751584"),])
unique(manifest450k$CHR) # Build 37
unique(manifest450k$Chromosome_36) # Build 36

# These are SNPs, not CpGs. Names start with "rs"
rsRows <- manifest450k[manifest450k$Chromosome_36 == "",]
mean(rsRows$Chromosome_36 == rsRows$CHR) # both builds mark SNPs the same way

# If we remove the SNP rows, then manifest450k has the same nrows as Locations
# Let's compare the order of probes across the two annotation sources
# starting with names
sum(rownames(Locations) == manifest450k$IlmnID[manifest450k$Chromosome_36 != ""])
# not in the same order at all
# Let's sort them
sortedLocations <- Locations[order(rownames(Locations)),]
sortedManifest <- manifest450k[order(manifest450k$IlmnID),]

# The two builds are definitely different
mean(manifest450k$MAPINFO[manifest450k$Chromosome_36 != ""] == manifest450k$Coordinate_36[manifest450k$Chromosome_36 != ""])

# Which build matches Locations, if any?
# Need to get rid of the "chr" string in Locations chromosome names
sortedLocations$chrNum <- sub("chr", "", sortedLocations[, "chr"])

# Probably only mismatched with "MULTI" in Chromosome_36
mean(sortedLocations$chrNum == sortedManifest$Chromosome_36[sortedManifest$Chromosome_36 != ""])
View(sortedManifest[which(sortedLocations$chrNum != sortedManifest$Chromosome_36[sortedManifest$Chromosome_36 != ""])])
mean(sortedLocations$pos == sortedManifest$Coordinate_36[sortedManifest$Chromosome_36 != ""])

mean(sortedLocations$chrNum == sortedManifest$CHR[sortedManifest$CHR != ""]) # perfect
mean(sortedLocations$pos == sortedManifest$MAPINFO[sortedManifest$CHR != ""]) # perfect
mean(Islands.UCSC[sortedManifest$IlmnID[sortedManifest$Chromosome_36 != ""], "Relation_to_Island"] == sortedManifest$Relation_to_UCSC_CpG_Island[sortedManifest$Chromosome_36 != ""])

Islands.UCSC["cg00000108", "Relation_to_Island"]
sortedManifest$Relation_to_UCSC_CpG_Island[sortedManifest$Relation_to_UCSC_CpG_Island == ""] <- "OpenSea"
mean(Islands.UCSC[sortedManifest$IlmnID[sortedManifest$Chromosome_36 != ""], "Relation_to_Island"] == sortedManifest$Relation_to_UCSC_CpG_Island[sortedManifest$Chromosome_36 != ""])

# Conclusion: use the build 37 columns, i.e. the CHR and MAPINFO columns
# CHR == Locations
# Relation_to_UCSC_CpG_Island == Islands.UCSC[,"Relation_to_Island"]
# MAPINFO == pos
neededColumns450k <- c("IlmnID", "CHR", "Relation_to_UCSC_CpG_Island", "MAPINFO", 
                       "Probe_SNPs_10", "Probe_SNPs")

essentialManifest450k <- manifest450k[, ..neededColumns450k]
# fwrite(essentialManifest450k, file = "/home/lutiffan/github/methylModes/IlluminaManifest450k.csv")
test450k <- fread(file = "/nfs/turbo/bakulski1/People/lutiffan/github/methylModes/IlluminaManifest450k.csv")

manifestEpic <- fread("/nfs/turbo/bakulski1/People/lutiffan/dataCleaningIlluminaManifests/MethylationEPIC_v-1-0_B4.csv", skip = 7)
neededColumnsEpic <- c("IlmnID", "CHR", "Relation_to_UCSC_CpG_Island", "MAPINFO", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")
essentialManifestEpic <- manifestEpic[1:865918, ..neededColumnsEpic]

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(manifestEpic, 20)
Locations[c("cg07881041", "cg18478105", "cg01763666"),] # Spot check
Islands.UCSC[c("cg07881041", "cg18478105", "cg01763666"), "Relation_to_Island"]
# EPIC annotation library matches build 37
# fwrite(essentialManifestEpic, file = "/home/lutiffan/github/methylModes/IlluminaManifestEPIC.csv")
testEpic <- fread(file = "/home/lutiffan/github/methylModes/IlluminaManifestEPIC.csv")

# Note that IlmnID is identical to Name for 450k, but not always so for EPIC
# In EPIC, the probes for which IlmnID != Name are not "cg..." probes
which(!(manifestEpic$IlmnID %like% "cg") & (manifestEpic$IlmnID != manifestEpic$Name))
manifestEpic$IlmnID[865919]; manifestEpic$Name[865919] # Exclude this row and beyond. Those are controls
sum(manifestEpic$IlmnID != manifestEpic$Name) # 636
length(which(!(manifestEpic$IlmnID %like% "cg") & (manifestEpic$IlmnID != manifestEpic$Name))) # 636