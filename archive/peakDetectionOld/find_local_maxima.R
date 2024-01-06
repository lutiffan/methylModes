# https://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset
# Problem: this code assumes that the slope is never flat

argmax <- function(x, y, w=1, ...) {
  require(zoo)
  #browser()
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, 
                     align="center")
  delta <- y.max - y.smooth[-c(1:w, (n+1)-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

test <- function(w, span) {
  peaks <- argmax(x, y, w=w, span=span)
  
  plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", 
              span = ", span, sep=""))
  lines(x, peaks$y.hat,  lwd=2) #$
  y.min <- min(y)
  sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, 
                                                    peaks$y.hat[i]),
                                    col="Red", lty=2))
  points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, 
         cex=1.25)
}

x <- 1:1000 / 100 - 5
y <- exp(abs(x)/20) * sin(2 * x + (x/5)^2) + cos(10*x) / 5 + 
  rnorm(length(x), sd=0.05)
par(mfrow=c(3,1))
test(2, 0.05)
test(30, 0.05)
test(2, 0.2)


# Run this after running above code
dev.off()

# https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
# Problem: struggles with small perturbations in flat parts of graph
# Method with diff()
tt <- c(1,2,3,2,1, 1, 2, 1)
which(diff(sign(diff(tt)))==-2)+1


localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  # rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

# Mock up some small examples of code that looks like of like beta distribtions
# Unimodal
testdat1 <- rbeta(1000, 2,2)
testdensity1 <- density(testdat1, from = 0, to = 1, adjust = 1.5)
plot(testdensity1$x, testdensity1$y)
which.max(testdensity1$y) # The easy answer
localMaxima(testdensity1$y) # works ok
which(diff(sign(diff(testdensity1$y)))==-2)+1 # works ok

# Bimodal
set.seed(10)
testdat2 <- rbeta(1000, 0.8,0.8)
testdensity2 <- density(testdat2, from = 0, to = 1, adjust = 1.5)
plot(testdensity2$x, testdensity2$y)
localMaxima(testdensity2$y) # works ok

# Trimodal
set.seed(7)
testdat3 <- rbeta(1000, 0.8,0.8)
testdensity3 <- density(testdat3, from = 0, to = 1, adjust = 1.5)
plot(testdensity3$x, testdensity3$y)
localMaxima(testdensity3$y) # works ok
which(diff(sign(diff(testdensity3$y)))==-2)+1 # works ok

testdat <- rbeta(1000, 2, 100)
testdensity <- density(testdat, from = 0, to = 1, adjust = 1.5)
# Try rounding to get rid of small values close to zero
testdensity$y <- ifelse(testdensity$y < 0.001, 0, testdensity$y)
plot(testdensity$x, testdensity$y)
which.max(testdensity$y) # The one and only right answer
argmax(x = testdensity$x, y = testdensity$y) # doesn't work well
localMaxima(testdensity$y) # also fails, but does better than argmax()
which(diff(sign(diff(testdensity$y)))==-2)+1 # same as localMaxima

smoothed <- loess(testdensity$y ~ testdensity$x)
plot(x = testdensity$x, smoothed$fitted)

# I don't know how to quickly simulate this kind of distribution
nonzero <- testdensity$y > 0.001

zeroRun1 <- seq(from = max(testdensity$x[nonzero] + 0.01), 
                to = max(testdensity$x[nonzero]) + 0.20, by = 0.01)
offset1 <- max(zeroRun1) + 0.01
compress1 <- 0.3 # Try making the peak shorter
zeroRun2 <- seq(from = max(testdensity$x[nonzero] + offset1) + 0.01,
                to = max(testdensity$x[nonzero] + offset1) + 0.20, 
                by = 0.01)
offset2 <- max(zeroRun2) + 0.01
compress2 <- 0.01
faketrimodalx <- c(testdensity$x[nonzero], 
                   zeroRun1,
                   testdensity$x[nonzero] + offset1, 
                   zeroRun2,
                   testdensity$x[nonzero] + offset2
                   )
faketrimodaly <- c(testdensity$y[nonzero], rep(0, 20), testdensity$y[nonzero]*compress1, 
                   rep(0, 20), testdensity$y[nonzero]*compress2)
plot(faketrimodalx, faketrimodaly)
which(diff(sign(diff(faketrimodaly)))==-2)+1
