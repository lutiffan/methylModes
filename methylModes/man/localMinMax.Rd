\name{localMinMax}
\alias{localMinMax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
localMinMax(fitted, zeroThreshold = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitted}{
%%     ~~Describe \code{fitted} here~~
}
  \item{zeroThreshold}{
%%     ~~Describe \code{zeroThreshold} here~~
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
function (fitted, zeroThreshold = 1e-06) 
{
    fitted <- ifelse(fitted < zeroThreshold, 0, fitted)
    plusMinus <- diff(c(-Inf, fitted)) > 0L
    slopeZero <- cumsum(rle(plusMinus)$lengths)
    slopeZero
    if (length(slopeZero) == 1) {
        return(list(maximaIdx = slopeZero, minimaIdx = c(1, 500)))
    }
    whereMax <- seq.int(1L, length(slopeZero), 2L)
    whereMin <- (1:length(slopeZero))[-whereMax]
    if (fitted[[1]] == fitted[[2]]) {
        whereMax <- whereMax[-1]
        whereMin <- whereMin[-1]
    }
    leftMin <- whereMax - 1
    leftMin <- leftMin[-1]
    rightMin <- whereMax + 1
    if (rightMin[length(rightMin)] > length(slopeZero)) {
        rightMin[length(rightMin)] = length(slopeZero)
    }
    data.frame(maximaIdx = slopeZero[whereMax], leftMinIdx = c(1, 
        slopeZero[leftMin]), rightMinIdx = slopeZero[rightMin])
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
