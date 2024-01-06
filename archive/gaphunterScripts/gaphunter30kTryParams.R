library(minfi)
library(gridExtra)
library(ggplot2)

# Read beta matrix with cross-reactive probes, teenagers, and sex chromosome probes removed
betas <- readRDS("cleanedBetasChildOnlyAutosomal.RDS")
# cleansnpbeta <- readRDS("cleansnpbeta.RDS")

thresholds <- c(0.025, 0.05, 0.10, 0.2)
outCutoffs <- c(0.005, 0.01, 0.05, 0.1)

allPlots <- vector(mode = "list", length = length(thresholds)*length(outCutoffs))

for (x in seq_along(thresholds)) {
  print(x)
  for (y in seq_along(outCutoffs)) {
    print(y)
    starttime <- Sys.time()
    ghresults <- gaphunter(object = betas, threshold = thresholds[x], 
                           outCutoff = outCutoffs[y])
    Sys.time() - starttime
    saveRDS(ghresults, file = paste0("/home/lutiffan/sensitivityAnalysis/gaphunter_threshold_", 
                                     gsub('\\.', '_', thresholds[x]),
                                     "_outCutoff_",
                                     gsub('\\.', '_', outCutoffs[y]), ".RDS"))
    allPlots[[4*(x - 1) + y]] <- ggplot(data = ghresults$proberesults, aes(Groups)) +
      geom_bar() +
      geom_text(stat='count', aes(label=..count..), vjust = -1) + 
      ggtitle(label = paste("t =", thresholds[x], "oC=", outCutoffs[y]))
  }
}

png(filename = "/home/lutiffan/googleSlides/gaphunterBetasBarplot.png", width = 600, height = 500)
grid.arrange(
  allPlots[[1]], allPlots[[2]], allPlots[[3]], allPlots[[4]], 
  allPlots[[5]], allPlots[[6]], allPlots[[7]], allPlots[[8]], 
  allPlots[[9]], allPlots[[10]], allPlots[[11]], allPlots[[12]],
  allPlots[[13]], allPlots[[14]], allPlots[[15]], allPlots[[16]],
  nrow = 4
)
dev.off()
