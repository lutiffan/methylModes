library(mixtools)

# Example data
realmu <- c(0, 4, 8)
realsd <- c(1, 1.5, 0.5)
data <- rnorm(1000, mean = realmu, sd = realsd)

# Fit Gaussian mixture model
fit <- normalmixEM(data, k = 1)

# Print the estimated parameters
print(fit)

# Plot the data and fitted mixture model
plot(fit)

mean1 <- 0.05
mean2 <- 0.1
mean3 <- 0.3
mean4 <- 0.6

sd1 <- 0.008
sd2 <- 0.03
sd3 <- 0.02
sd4 <- 0.02

set.seed(30)
peak1 <- rnorm(900, mean =  mean1, sd = sd1)
peak2 <- rnorm(200, mean = mean2, sd = sd2)
peak3 <- rnorm(50, mean = mean3, sd = sd3)
peak4 <- rnorm(500, mean = mean4, sd = sd4)
fourPeaks <- c(peak1, peak2, peak3, peak4)

fit2 <- normalmixEM(fourPeaks, k = 4)
fit2$mu; fit2$sigma


####  Iterative  ####

library(mixtools)
library(MASS)

# Example data
randomdata <- rnorm(1000, mean = c(0, 4, 8), sd = c(1, 1.5, 0.5))

# Constants
MAX_K <- 4
LL_INC_THRESHOLD <- 50
# BIC_INC_THRESHOLD <- 100

# Initialize variables
ki <- 1
BEST_MODEL <- NULL
#BEST_BIC <- Inf
BEST_LL <- Inf

while (ki <= MAX_K) {
  # Fit Gaussian mixture model with k components
  if (ki == 1) {
    fit <- fitdistr(fourPeaks, "normal")
    # fit <- normalmixEM(x = randomdata, lambda = NULL, mu = NULL, sigma = NULL, 
    #                    k = 1, 
    #             mean.constr = NULL, sd.constr = NULL,
    #             epsilon = 1e-08, maxit = 1000, maxrestarts = 20, 
    #             verb = FALSE, fast = FALSE, ECM = FALSE,
    #             arbmean = TRUE, arbvar = TRUE)
    # 
    # # Extracted from inside normalmixEM
    # normalmix.init(x = randomdata, lambda = NULL, mu = NULL, s = NULL, 
    #                 k = 1, arbmean = FALSE, arbvar = TRUE)
  } else {
    fit <- normalmixEM(x = randomdata, k = ki)
  }
  
  # Calculate BIC for the current model
  # bic <- BIC(fit)
  ll <- fit$loglik
  
  if (ki == 1 || ll > BEST_LL + LL_INC_THRESHOLD) {
    BEST_MODEL <- fit
    BEST_LL <- ll
  } else {
    break
  }
  
  print(ll)
  # if (ki == 1 || bic < BEST_BIC + BIC_INC_THRESHOLD) {
  #   BEST_MODEL <- fit
  #   BEST_BIC <- bic
  # } else {
  #   break
  # }
  
  ki <- ki + 1
}

# Print the best model and its BIC
cat("Best Model (k =", length(BEST_MODEL$lambda), "):", "\n")
print(BEST_MODEL)
cat("BIC:", BEST_BIC, "\n")

# Plot the data and best model
plot(data, density = 30)
lines(density(data))
lines(density(rmixtools:::densityfun(BEST_MODEL, "x")), col = "red")

regmixmodel.sel(x = randomdata)







library(mixtools)

# Example data
data <- rnorm(1000, mean = c(0, 4, 8), sd = c(1, 1.5, 0.5))
dataPts <- rnorm(1000, mean = 0.5, sd = 0.1)

# Constants
MAX_K <- 2
BIC_INC_THRESHOLD <- 100

# Initialize variables
k <- 1
BEST_MODEL <- NULL
BEST_BIC <- Inf

while (k <= MAX_K) {
  # Fit Gaussian mixture model with k components
  if (k == 1) {
    fit <- fitdistr(dataPts, "normal")
    bic <- BIC(fit)
  } else {
    fit <- normalmixEM(dataPts, k = k)
    # Extract log-likelihood and number of parameters
    loglik <- fit$loglik
    n <- length(data)
    p <- k * 3 - 1  # Number of parameters: mean, standard deviation, and mixing proportion for each component
    
    # Calculate BIC for the current model
    bic <- -2 * loglik + p * log(n)
  }
  
  
  if (k == 1 || bic < BEST_BIC - BIC_INC_THRESHOLD) {
    BEST_MODEL <- fit
    BEST_BIC <- bic
  } else {
    break
  }
  
  k <- k + 1
}

BEST_MODEL$mu

# Print the best model and its BIC
cat("Best Model (k =", length(BEST_MODEL$lambda), "):", "\n")
print(BEST_MODEL)
cat("BIC:", BEST_BIC, "\n")

# Plot the data and best model
plot(data, density = 30)
lines(density(data))
lines(density(rmixtools:::densityfun(BEST_MODEL, "x")), col = "red")
