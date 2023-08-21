idats <- list.files("/nfs/turbo/bakulski1/People/johndou/Fragile_Families/IDATs/", 
                    recursive = T, full.names = T, pattern = ".idat$")
idats <- idats %>% gsub("_Grn.idat|_Red.idat", "", .) %>% .[duplicated(.)]
meth <- read_idats(idats) %>% detectionP 
# saveRDS(meth, "final/pp/meth_withdetp.rds")

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
pd <- data.table(pd,ctrls)

#Sex check
pd[, c("x","y") := meth %>% correct_dye_bias %>% check_sex]
pd[, predsex := predict_sex(x,y,which(sex == "m"), which(sex == "f"))]

#Create beta matrix for SNP probes and use to create genotype data
snps <- meth$manifest$probe_id[startsWith(meth$manifest$probe_id, "rs")]
toolsbeta <- meth %>% mask(0.01) %>% correct_dye_bias() %>% dont_normalize
snpbeta <- toolsbeta[rownames(toolsbeta) %in% snps,]

geno <- call_genotypes(snpbeta, learn = T) #Learn=T sets EM algorithm to adapt to our data
# saveRDS(geno, "final/pp/called_genotypes.rds")
pd[, snpout := snp_outliers(geno)]
agree <- check_snp_agreement(geno, pd[, idnum], pd[, MethID])

rm(meth)

failedprobes <- readRDS("final/pp/failed_probes.rds")
badprobes <- names(failedprobes)[failedprobes > 0.05]
str(badprobes)