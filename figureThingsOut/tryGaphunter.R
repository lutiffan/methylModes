# GHresults30k <- gaphunter(object = betas[1:30000,], verbose = TRUE)
library(ggplot2)
library(VennDiagram)
library(grid)
source("visualizeBetaPeaks.R")
GHresults30k <- readRDS("GHresults30k_defaultParams.RDS")
# 1020 gap probes found

peakSummary0105 <- readRDS("peakSummary_like_gaphunter.RDS")
# The object "peakSummary" is hard-coded into visualization functions
# peakSummary <- peakSummary0105
# proportionSample <- 0.01
# personalSpace <- 0.05

peaksummary0505 <- readRDS("sensitivityAnalysis/peakSummary_sensitivity_prop_0_05_space_0_05_start_1_end_30000.RDS")
# Parameters associated with this peakSummary object
kernelType = "gaussian" # "This must partially match one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine", with default "gaussian", and may be abbreviated to a unique prefix (single letter)."
bandwidthType = NA # options: NA, "sheatherJones"
numBreaks = 500 # The "n" parameter for density(): "the number of equally spaced points at which the density is to be estimated..." You need at least 500 for a decent smooth curve
densityAdjust = 1.5
pushToZero = 0.0001 # Get rid of minuscule floating point numbers  created by density() that really should be zero

rangeStart = 1
rangeEnd = 30000

betaFileName <- "/nfs/turbo/bakulski1/People/data/ffcw/methyl_data/Feb2022_datafreeze/beta_analytic_freeze1.rds" # file.choose()
betas <- readRDS(betaFileName) # 5.9 GB

numGroupsSum <- table(GHresults30k$proberesults$Groups); numGroupsSum
barplot(numGroupsSum)

# heresadataframe <- data.frame("nGroups" = nGroups)
# ggplot(data = heresadataframe, aes(groups)) +
#   geom_bar() +
#   geom_text(stat='count', aes(label=..count..), vjust = -1)

# Examine intersection with my algorithm
ghprobes <- rownames(GHresults30k$proberesults)

##### Comparison 1: all "gap" probes #####
multiIdx <- which(peakSummary$numPeaks >= 2) + rangeStart - 1
myprobes <- rownames(betas)[multiIdx]
overlap <- intersect(myprobes, ghprobes)
# Which gaphunter probes do I miss?
ghOnly <- setdiff(ghprobes, myprobes)
mineOnly <- setdiff(myprobes, ghprobes)

for (m in seq_along(ghOnly)) {
  peakSummaryPlot(probe.id = ghOnly[m], range.start = rangeStart,
                  range.end = rangeEnd)
  ghcount <- GHresults30k$proberesults$Groups[which(ghprobes == ghOnly[m])]
  legend("top",     as.character(ghcount), pch=1, title= "Gaphunter count: ", inset = .05)
  browser()
}

for (m in seq_along(mineOnly)) {
  peakSummaryPlot(probe.id = mineOnly[m], range.start = rangeStart,
                  range.end = rangeEnd)
  ghcount <- GHresults30k$proberesults$Groups[which(ghprobes == mineOnly[m])]
  legend("top",     as.character(ghcount), pch=1, title= "Gaphunter count: ", inset = .05)
  browser()
}

# Mine doesn't look quite right:
# cg00995284
# peakSummaryPlot(probe.id = "cg00995284", range.start = rangeStart, range.end = rangeEnd)
seePeaks(probe.id = "cg00995284")
# cg16348638
seePeaks(probe.id = "cg16348638")
which(rownames(betas) == "cg16348638") # 20284
# cg20268341
seePeaks(probe.id = "cg20268341")
which(rownames(betas) == "cg20268341") # 28287

# big contrast:
# cg27278382
# cg22403382

# trimodal????!!
# cg25950377

venn.diagram(list("GH" = ghprobes, "myProbes" = myprobes), 
             category.names = c("Gaphunter", "Sophia"),
             filename = paste0("googleSlides/vennDiagrams/", "gaphunter_venn_", 
                               "allgroups",
                               "prop_", proportionSample, 
                               "_space_", personalSpace),
             disable.logging = TRUE,
             main = "All gap probes identified by Gaphunter and Sophia's algorithm",
             width = 3500)

##### Comparison 2: bi-modal probes #####
bimodalIdx <- which(peakSummary$numPeaks == 2) + rangeStart - 1
myProbes2 <- rownames(betas)[bimodalIdx]
ghProbes2 <- rownames(GHresults30k$proberesults[GHresults30k$proberesults$Groups == 2,])
venn.diagram(list("GH2" = ghProbes2, "myProbes2" = myProbes2), 
             category.names = c("Gaphunter", "Sophia"),
             filename = paste0("googleSlides/vennDiagrams/", "gaphunter_venn_",
                               "2groups",
                               "prop_", proportionSample, 
                               "_space_", personalSpace),
             disable.logging = TRUE,
             main = "Two-group gap probes identified by Gaphunter and Sophia's algorithm",
             width = 3500)

overlap2GroupCounts <- GHresults30k$proberesults$Groups[rownames(GHresults30k$proberesults) %in% myProbes2]
table(overlap2GroupCounts)

ghOnly2 <- setdiff(ghProbes2, myProbes2)
# set.seed(7323)
# biModalSample <- sample(1:length(ghOnly2), 60, replace = F)
for (m in seq_along(ghOnly2)) {
  peakSummaryPlot(probe.id = ghOnly2[m], range.start = rangeStart,
                  range.end = rangeEnd)
  ghcount <- GHresults30k$proberesults$Groups[which(ghprobes == ghOnly2[m])]
  legend("top",     as.character(ghcount), pch=1, title= "Gaphunter count: ", inset = .05)
  browser()
}

# Mine doesn't look quite right
# cg00997411
seePeaks(probe.id = "cg00997411")
# cg18107019
seePeaks(probe.id = "cg18107019")
# cg00995284
seePeaks(probe.id = "cg00995284")
# cg25320763

# Ok that's a subtle difference
# cg27558057
# cg06176944

##### Comparison 3: tri-modal probes #####
trimodalIdx <- which(peakSummary$numPeaks == 3) + rangeStart - 1
myProbes3 <- rownames(betas)[trimodalIdx]
ghProbes3 <- rownames(GHresults30k$proberesults[GHresults30k$proberesults$Groups == 3,])
length(myProbes3)
length(ghProbes3)
venn.diagram(list("GH3" = ghProbes3, "myProbes3" = myProbes3), 
             category.names = c("Gaphunter", "Sophia"),
             filename = paste0("googleSlides/vennDiagrams/", "gaphunter_venn_",
                               "3groups",
                               "prop_", proportionSample, 
                               "_space_", personalSpace),
             disable.logging = TRUE,
             main = "Three-group gap probes identified by Gaphunter and Sophia's algorithm",
             width = 3500)



##### 
manymodalIdx <- which(peakSummary$numPeaks > 3) + rangeStart - 1
myProbes4 <- rownames(betas)[manymodalIdx]
ghProbes4 <- rownames(GHresults30k$proberesults[GHresults30k$proberesults$Groups > 3,])
length(myProbes4)
length(ghProbes4)
venn.diagram(list("GH4" = ghProbes4, "myProbes4" = myProbes4), 
             category.names = c("Gaphunter", "Sophia"),
             filename = paste0("googleSlides/vennDiagrams/", "gaphunter_venn_",
                               "manygroups",
                               "prop_", proportionSample, 
                               "_space_", personalSpace),
             disable.logging = TRUE,
             main = "Over three-group gap probes identified by Gaphunter and Sophia's algorithm",
             width = 3500)

countDiff <- data.frame("probe" = overlap, 
                        "Gaphunter" = numeric(length(overlap)), 
                        "mine" = numeric(length(overlap)))
for (i in seq_along(overlap)) {
  peakSummaryPlot(probe.id = overlap[i], range.start = rangeStart, 
                  range.end = rangeEnd)
 
  ghcount <- GHresults30k$proberesults$Groups[which(ghprobes == overlap[i])]
  mycount <- peakSummary$numPeaks[which(rownames(betas) == overlap[i])]
  
  countDiff[i,"GH"] <- ghcount
  countDiff[i, "mine"] <- mycount
  
  # if (ghcount != mycount) {
  #   print(overlap[i])
  #   print(paste("gaphunter:", GHresults30k$proberesults$Groups[which(ghprobes == overlap[i])]))
  #   print(paste("sophia:", peakSummary$numPeaks[which(rownames(betas) == overlap[i])]))
  # }
}
sum(countDiff$mine != countDiff$GH)
minelessgh <- countDiff[countDiff$mine < countDiff$GH,]; minelessgh
peakSummaryPlot(probe.id = "cg07136920", range.start = rangeStart, 
                range.end = rangeEnd)
peakSummaryPlot(probe.id = "cg19961153", range.start = rangeStart, 
                range.end = rangeEnd)

minemoregh <- countDiff[countDiff$mine > countDiff$GH,]; minemoregh

for (j in 1:nrow(minemoregh)) {
  peakSummaryPlot(probe.id = minemoregh[j,"probe"], range.start = rangeStart, 
                  range.end = rangeEnd)
  browser()
}

# GH probes that I don't catch
ghOnly <- setdiff(ghprobes, myprobes)
for (g in 1:length(ghOnly)) {
  peakSummaryPlot(probe.id = ghOnly[g], range.start = rangeStart, 
                  range.end = rangeEnd)
  ghcount <- GHresults30k$proberesults$Groups[which(ghprobes == ghOnly[g])]
  print(paste("gaphunter count:", ghcount))
  mycount <- peakSummary$numPeaks[which(rownames(betas) == ghOnly[g])]
  print(paste("my count:", mycount))
  browser()
}

# interesting results:
# cg25320763 doesn't look right (row 23139)

# barplot(height = numGroupsSum)

# saveRDS(GHresults30k, file = "GHresults30k_defaultParams.RDS")

# Crashes with 15 GB memory
# GHall <- gaphunter(object = betas, verbose = TRUE)
