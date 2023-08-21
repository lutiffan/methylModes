length(checkCrowding)
rle(checkCrowding)$lengths + 1

# diff vector begins with F
test <- c(F, F, T, T ,T ,F, F, T, T ,F); test
which(test)
which(diff(test) == 1)
clustered <- sort(c(which(test), which(diff(test) == 1))) + 1; clustered
testHeights <- c(0,0,1,2,3,2,0,1,2,3,0)

# diff vector begins with T
test2 <- c(T, T, F, F, T, T, T)
which(test2)
testChanges2 <- diff(test2); testChanges2
c(1, sort(c(which(test2), which(diff(test2) == 1))) + 1 )

test <- c(T)
c(1, sort(c(which(test), which(diff(test) == 1))) + 1 )

# Mark all crowded peaks as "not the tallest in its cluster" by default
clusterTallest <- logical(length(clustered))
clusterTallestIdx <- 1
clusterTallestHt <- testHeights[clustered[1]]; clusterTallestHt
numCrowded <- length(clustered)
# Keep track of which fake peak is at the start and end of clusters
numClusters <- sum(diff(clustered) >= 2) + 1
bookEnds <- data.frame("start" = integer(numClusters), "end" = integer(numClusters))
bookEnds$start[1] <- clustered[1]
clusterIdx <- 1
# Among consecutive crowded peaks, find the tallest
for (i in 2:numCrowded) {
  if ((clustered[i] - clustered[i - 1]) > 1) { 
    # Mark the end of the old cluster
    bookEnds$end[clusterIdx] <- clustered[i - 1]
    
    # Mark a cluster's tallest 
    clusterTallest[clusterTallestIdx] <- TRUE
    
    clusterTallestIdx <- clustered[i]
    clusterTallestHt <- testHeights[clustered[i]]
    
    clusterIdx <- clusterIdx + 1
    
    # Keep track of the start of the new cluster
    bookEnds$start[clusterIdx] <- clustered[i]
  }
  
  if (testHeights[clustered[i]] > clusterTallestHt) {
    clusterTallestIdx <- i
    clusterTallestHt <- testHeights[clustered[i]]
    
    if (i == numCrowded) {
      # Mark the end of the last cluster
      bookEnds$end[clusterIdx] <- clustered[i]
      
      clusterTallest[clusterTallestIdx] <- TRUE
    }
  }
  
  clusterTallestIdx; clusterTallestHt
}
clusterTallest