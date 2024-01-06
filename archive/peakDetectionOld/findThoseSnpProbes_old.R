# Needs 75+ GiB memory to load in meth object
# Need 130+ GB memory to run QC steps
library(dplyr)
library(minfi)
source("visualizeBetaPeaks.R")

snpbeta <- readRDS("~/snpbeta.RDS")

annotations <- read.csv(file = "/nfs/turbo/bakulski1/People/data/ffcw/methyl_data/Feb2022_datafreeze/pd_analytic_freeze1.csv") 
cleansnpbeta <- snpbeta[,colnames(snpbeta) %in% annotations$MethID]

numNAs <- numeric(nrow(cleansnpbeta))
for ( i in 1:nrow(cleansnpbeta)) {
  numNAs[i] <- sum(is.na(cleansnpbeta[i,]))
}

# Why do we still have NAs? Cut them for now
cleansnpbeta <- cleansnpbeta[numNAs == 0,]

ghSNPs <- gaphunter(object = cleansnpbeta, verbose = T)

ghCounts <- table(ghSNPs$proberesults$Groups)
barplot(ghCounts, main = "Gaphunter Group Counts Among 44 SNP Probes")

betas <- cleansnpbeta
for (i in 1:nrow(cleansnpbeta)) {
  seePeaks(row.id = i)
  legend("top", as.character(ghSNPs$proberesults$Groups[i]), pch=1, title= "Gaphunter group count: ", inset = .05)
  browser()
}

##### Ignore these QC steps below #####
library(ewastools)
library(data.table)

idats <- list.files("/nfs/turbo/bakulski1/People/johndou/Fragile_Families/IDATs/", 
                    recursive = T, full.names = T, pattern = ".idat$")
idats <- idats %>% gsub("_Grn.idat|_Red.idat", "", .) %>% .[duplicated(.)]
meth <- read_idats(idats) %>% detectionP 
# saveRDS(meth, "meth_withdetp.rds")

failedobs <- meth$detP > 0.01 | meth$N < 4 | meth$V < 4
failedobs[is.na(failedobs)] <- T # NA values assigned to failed observation

colnames(failedobs) <- meth$meta$sample_id
rownames(failedobs) <- meth$manifest$probe_id
# saveRDS(failedobs, "final/pp/failedobs.rds")

failedprobes <- rowMeans(failedobs)
failedsamps <- colMeans(failedobs[failedprobes < 0.05,]) # Pull out failed probes before calculating failed samps

# saveRDS(failedprobes, "final/pp/failed_probes.rds")
# saveRDS(failedsamps, "final/pp/failed_samps.rds")

#Controls
ctrls <- meth %>% control_metrics
# pd <- data.table(pd,ctrls)
pd <- data.table(ctrls)

saveRDS(pd, "pd.RDS")

#Sex check
newValues <- meth %>% correct_dye_bias %>% check_sex
pd[, c("x","y") := meth %>% correct_dye_bias %>% check_sex]
pd[, predsex := predict_sex(x,y,which(sex == "m"), which(sex == "f"))]

#Create beta matrix for SNP probes and use to create genotype data
snps <- meth$manifest$probe_id[startsWith(meth$manifest$probe_id, "rs")]
toolsbeta <- meth %>% mask(0.01) %>% correct_dye_bias() %>% dont_normalize
snpbeta <- toolsbeta[rownames(toolsbeta) %in% snps,]

# I haven't run this section yet
geno <- call_genotypes(snpbeta, learn = T) #Learn=T sets EM algorithm to adapt to our data
# saveRDS(geno, "final/pp/called_genotypes.rds")
pd[, snpout := snp_outliers(geno)]
agree <- check_snp_agreement(geno, pd[, idnum], pd[, MethID])

rm(meth)

failedprobes <- readRDS("final/pp/failed_probes.rds")
badprobes <- names(failedprobes)[failedprobes > 0.05]
str(badprobes)