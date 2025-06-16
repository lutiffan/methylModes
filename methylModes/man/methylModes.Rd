\name{methylModes}
\alias{methylModes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
methylModes(row.data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{row.data}{
%%     ~~Describe \code{row.data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (row.data = NULL) 
{
    if (is.na(bandwidthType)) {
        bandwidth = "nrd0"
    }
    else if (bandwidthType == "sheatherJones") {
        bandwidth = bw.SJ(row.data)
    }
    probeDensityEst <- density(row.data, from = 0, to = 1, n = numBreaks, 
        adjust = densityAdjust, bw = bandwidth, kernel = kernelType)
    detected <- localMinMax(probeDensityEst$y, zeroThreshold = pushToZero)
    nearCutoffPropSample <- NA
    nearCutoffPeakDistance <- NA
    if (nrow(detected) > 1) {
        leftMinx <- probeDensityEst$x[detected$leftMinIdx]
        rightMinx <- probeDensityEst$x[detected$rightMinIdx]
        propSample <- numeric(nrow(detected))
        for (p in seq_along(propSample)) {
            propSample[p] <- mean(row.data > leftMinx[p] & row.data < 
                rightMinx[p])
        }
        detected$propSample <- propSample
        belowCutoff <- propSample < proportionSample
        nearCutoffPropSample <- any(abs(propSample - proportionSample) <= 
            0.01)
        oneMinimaZero <- (probeDensityEst$y[detected$leftMinIdx] < 
            pushToZero) | (probeDensityEst$y[detected$rightMinIdx] < 
            pushToZero)
        detected <- detected[!(belowCutoff & oneMinimaZero), 
            ]
    }
    else {
        detected$propSample <- 1
    }
    xValues <- probeDensityEst$x[detected$maximaIdx]
    spacings <- diff(xValues)
    checkCrowding <- spacings <= peakDistance
    if (sum(checkCrowding) > 0) {
        nearCutoffPeakDistance <- any(abs(spacings - peakDistance) <= 
            0.01)
        clusters <- rle(checkCrowding)$lengths + 1
        if (checkCrowding[1]) {
            clusteredPeakIdx <- c(1, sort(c(which(checkCrowding), 
                which(diff(checkCrowding) == 1))) + 1)
            clusterSizes <- clusters[seq(from = 1, to = length(clusters), 
                by = 2)]
        }
        else {
            clusteredPeakIdx <- sort(c(which(checkCrowding), 
                which(diff(checkCrowding) == 1))) + 1
            clusteredPeakIdx
            clusterSizes <- clusters[seq(from = 2, to = length(clusters), 
                by = 2)]
        }
        clusterStart <- cumsum(clusterSizes) - clusterSizes + 
            1
        clusterEnd <- cumsum(clusterSizes)
        heights <- probeDensityEst$y[detected$maximaIdx[clusteredPeakIdx]]
        clusterTallest <- logical(length(clusteredPeakIdx))
        for (x in seq_along(clusterSizes)) {
            clusterWinner <- which.max(heights[clusterStart[x]:clusterEnd[x]]) + 
                clusterStart[x] - 1
            clusterTallest[clusterWinner] <- TRUE
            detected$leftMinIdx[clusteredPeakIdx[clusterWinner]] <- detected$leftMinIdx[clusteredPeakIdx[clusterStart[x]]]
            detected$rightMinIdx[clusteredPeakIdx[clusterWinner]] <- detected$rightMinIdx[clusteredPeakIdx[clusterEnd[x]]]
        }
        clearWinners <- setdiff(1:nrow(detected), clusteredPeakIdx)
        detected <- detected[c(clearWinners, clusteredPeakIdx[clusterTallest]), 
            ]
        detected <- detected[order(detected$maximaIdx), ]
    }
    if (length(detected$maximaIdx) > 1) {
        leftMinx <- probeDensityEst$x[detected$leftMinIdx]
        rightMinx <- probeDensityEst$x[detected$rightMinIdx]
        propSample <- numeric(nrow(detected))
        for (p in seq_along(propSample)) {
            propSample[p] <- mean(row.data > leftMinx[p] & row.data < 
                rightMinx[p])
        }
        detected$propSample <- propSample
        aboveCutoff <- propSample > proportionSample
        detected <- detected[aboveCutoff, ]
    }
    else {
        detected$propSample <- 1
    }
    mins <- unique(c(detected$leftMinIdx, detected$rightMinIdx))
    if (nrow(detected) > 1) {
        gap <- sum(probeDensityEst$y[mins[2:(length(mins) - 1)]] < 
            pushToZero) > 0
    }
    else {
        gap <- NA
    }
    list(detected = detected, probeDensityEst = probeDensityEst, 
        gap = gap, nearCutoffPropSample = nearCutoffPropSample, 
        nearCutoffPeakDistance = nearCutoffPeakDistance)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
