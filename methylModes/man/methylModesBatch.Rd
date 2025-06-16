\name{methylModesBatch}
\alias{methylModesBatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
methylModesBatch(betas = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betas}{
%%     ~~Describe \code{betas} here~~
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
function (betas = NULL) 
{
    if (is.null(betas)) 
        stop(simpleError("Invalid beta matrix."))
    template <- data.table(probeName = character(1), numPeaks = numeric(1), 
        meanBeta = numeric(1), peakLocations = vector(mode = "list", 
            length = 1), leftMin = vector(mode = "list", length = 1), 
        rightMin = vector(mode = "list", length = 1), proportionSample = vector(mode = "list", 
            length = 1), peakVariance = vector(mode = "list", 
            length = 1))
    cl <- parallelly::makeClusterPSOCK(availableCores(omit = 1), 
        autoStop = TRUE)
    doParallel::registerDoParallel(cl)
    peakSummary <- foreach::foreach(probe = iter(betas, by = "row"), 
        .combine = "rbind", .packages = c("foreach", "data.table"), 
        .export = c("methylModes", "localMinMax", "proportionSample", 
            "peakDistance", "kernelType", "bandwidthType", "numBreaks", 
            "densityAdjust", "pushToZero")) \%dopar\% {
        foundPeaks <- methylModes(row.data = probe)
        detected <- foundPeaks$detected
        probeDensityEst <- foundPeaks$probeDensityEst
        peakVariance = numeric(nrow(detected))
        for (k in 1:nrow(detected)) {
            leftBeta <- probeDensityEst$x[detected$leftMinIdx[k]]
            rightBeta <- probeDensityEst$x[detected$rightMinIdx[k]]
            peakVariance[k] <- var(probe[((probe > leftBeta) & 
                (probe < rightBeta))])
        }
        meanBeta <- mean(probe)
        template[1, ] <- list(probeName = rownames(probe), numPeaks = nrow(detected), 
            meanBeta = meanBeta, peakLocations = probeDensityEst$x[detected$maximaIdx], 
            leftMin = probeDensityEst$x[detected$leftMinIdx], 
            rightMin = probeDensityEst$x[detected$rightMinIdx], 
            proportionSample = detected$propSample, peakVariance = peakVariance)
        template
    }
    stopCluster(cl)
    peakSummary
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
