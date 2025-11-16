methylToSNPFixed<- function (data, gap.ratio = 0.75, gap.sum.ratio = 0.5, verbose = FALSE, 
          outlier.sd = 3, SNP) 
{
  if ((!length(dim(data)) == 2) || (dim(data)[2] < 2)) {
    stop("[MethylToSNP] There have to be at least 2 samples")
  }
  if (dim(data)[2] < 50) {
    message("[MethylToSNP] Warning, SNP detection may be unreliable in datasets with less than 50 samples")
  }
  if (gap.ratio >= 1 || gap.ratio <= 0) {
    stop("[MethylToSNP] You have entered an unacceptable gap.ratio value. This must be a number between zero and one.")
  }
  if (gap.sum.ratio <= 0 || gap.sum.ratio >= 1) {
    stop("[MethylToSNP] You have entered an unacceptable gap.sum.ratio value. This must be a number between zero and one.")
  }
  if (is(data, "GenomicRatioSet") || is(data, "GenomicMethylSet") || 
      is(data, "MethylSet") || is(data, "RatioSet")) {
    if (verbose) 
      message("[MethylToSNP] Calculating beta matrix. minfi is required")
    .checkMinfi()
    data <- getBeta(data)
  }
  if (missing(SNP)) {
    message("[MethylToSNP] Optionally, specify SNPs in a data frame with row names corresponding to cg probes (such as SNPs.147CommonSingle in minfiData or minfiDataEPIC package)")
  }
  potentialSNPs <- NULL
  probes <- rownames(data)
  for (i in 1:dim(data)[1]) {
    probe <- probes[i]
    x <- as.numeric(na.omit(data[i, ]))
    if (length(x) < 3) {
      next
    }
    span <- diff(range(x, na.rm = TRUE))
    if (span >= 0.5) {
      weights <- 1/approxfun(density(x))(x)
      kmeans <- Ckmeans.1d.dp(x = t(x), y = weights, k = 3)
      clusters <- kmeans$cluster
      centers <- kmeans$centers
      top <- which(centers == max(centers))
      bottom <- which(centers == min(centers))
      middle <- c(1:3)[-c(top, bottom)]
      if (outlier.sd != FALSE) {
        top_outliers <- which(abs(scale(x[clusters == 
                                            top])) > outlier.sd)
        clusters[top_outliers] <- 4
        middle_outliers <- which(abs(scale(x[clusters == 
                                               middle])) > outlier.sd)
        clusters[middle_outliers] <- 4
        bottom_outliers <- which(abs(scale(x[clusters == 
                                               bottom])) > outlier.sd)
        clusters[bottom_outliers] <- 4
        
        if (length(x[clusters == top]) == 0 | 
            length(x[clusters == middle]) == 0 | 
            length(x[clusters == bottom]) == 0) {
          next
        }
        
        # if (length(x[clusters == top]) * length(x[clusters == 
        #                                           middle]) * length(x[clusters == bottom]) == 
        #     0) {
        #   next
        # }
      }
      gaps.sorted <- sort(c(min(x[clusters == top]) - max(x[clusters == 
                                                              middle]), min(x[clusters == middle]) - max(x[clusters == 
                                                                                                             bottom])), decreasing = TRUE)
      gaps.largest <- gaps.sorted[1]
      gaps.smallest <- gaps.sorted[2]
      if ((gaps.smallest >= gap.ratio * gaps.largest) && 
          (sum(gaps.sorted) >= gap.sum.ratio * span)) {
        if (verbose) {
          message(probe)
        }
        potentialSNPs <- append(potentialSNPs, probe)
      }
    }
    else {
    }
    if (verbose) {
      if ((i%%1000) == 0) {
        message("[MethylToSNP] Processed: ", i, " Identified potential SNPs: ", 
                length(potentialSNPs))
      }
    }
  }
  if (verbose) {
    if (length(potentialSNPs) > 0) {
      message("[MethylToSNP] Number of potential SNPs found: ", 
              length(potentialSNPs))
    }
    else {
      warning("[MethylToSNP] No potential SNPs found.")
    }
  }
  snp.conf <- rep(NA, length(potentialSNPs))
  samples.low <- rep(NA, length(potentialSNPs))
  samples.mid <- rep(NA, length(potentialSNPs))
  samples.high <- rep(NA, length(potentialSNPs))
  for (i in 1:length(snp.conf)) {
    count.low <- length(which(data[potentialSNPs[i], ] <= 
                                0.25))
    count.high <- length(which(data[potentialSNPs[i], ] >= 
                                 0.75))
    count.mid <- length(data[potentialSNPs[i], ]) - count.low - 
      count.high
    L <- (count.high > 0) * (count.mid > 0) * (count.low > 
                                                 0)
    snp.conf[i] <- round(L * (count.high + (count.mid/2))/dim(data)[2], 
                         2)
    samples.low[i] <- count.low
    samples.mid[i] <- count.mid
    samples.high[i] <- count.high
  }
  if (length(potentialSNPs) > 0) {
    results <- data.frame(row.names = potentialSNPs, confidence = snp.conf, 
                          samples_low = samples.low, samples_mid = samples.mid, 
                          samples_high = samples.high)
    if (!missing(SNP)) {
      results <- cbind(results, SNP[potentialSNPs, ])
    }
    return(results)
  }
}
