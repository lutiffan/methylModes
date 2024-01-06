library(mclust)
source("visualizeBetaPeaks.R")
source("standardParams.R")

betas <- readRDS("cleanedBetasChildOnlyAutosomal.RDS")

#### GMMC ####

row1 <- betas[1,]
fit1 <- Mclust(data = row1)
G <- fit1$G
classes <- fit1$classification
classRanges <- vector(mode = "list", length = G)
for (i in 1:G) {
  inClass <- classes == i
  classRanges[[i]] <- range(row1[inClass])
}

color = c("orange", "purple", "cyan")
hist(row1, breaks = 500)
for (i in 1:G) {
  abline(v = classRanges[[i]][1], col = color[i], lwd = 2)
  abline(v = classRanges[[i]][2], col = color[i], lwd = 2)
}

seeGMMC <- function(data, groups, classification, 
                    colors = c("orange", "purple", "cyan", "darkred",
                               "bisque4", "darkgoldenrod3", "aquamarine4",
                               "coral2")) {
  classRanges <- vector(mode = "list", length = length(groups))
  
  for (i in seq_along(groups)) {
    inClass <- classification == groups[i]
    classRanges[[i]] <- range(data[inClass])
  }
  #browser()
  for (i in seq_along(groups)) {
    abline(v = classRanges[[i]][1], col = colors[i], lwd = 2)
    abline(v = classRanges[[i]][2], col = colors[i], lwd = 2)
  }
}

seeGH <- function(data, group.num, classification, 
                         colors = c("orange", "purple", "cyan", "darkred",
                                    "bisque4", "darkgoldenrod3", "aquamarine4",
                                    "coral2")) {
  classRanges <- vector(mode = "list", length = group.num)
  
  for (i in 1:group.num) {
    inClass <- classification == i
    classRanges[[i]] <- range(data[inClass])
  }
  #browser()
  for (i in 1:group.num) {
    abline(v = classRanges[[i]][1], col = colors[i], lwd = 2)
    abline(v = classRanges[[i]][2], col = colors[i], lwd = 2)
  }
}

seePeaks(row.id = 1)
# plot(fit1, what = "density")
seeGMMCandGH(data = row1, group.num = fit1$G, classification = fit1$classification)

#########

#### Gaphunter ####
betaIdx <- which(rownames(betas) %in% rownames(GHresults$proberesults))[1] # 51
ghIdx <- which(rownames(GHresults$sampleresults) == rownames(betas)[betaIdx])
ghGroups <- GHresults$proberesults$Groups[ghIdx] # Number of groups
ghClass <- GHresults$sampleresults[ghIdx,] # Classification for this probe

seePeaks(probe.id = rownames(betas)[betaIdx])
seeGMMCandGH(data = betas[betaIdx,], group.num = ghGroups, classification = ghClass)
