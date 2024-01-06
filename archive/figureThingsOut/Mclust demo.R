library(mclust)

# Example data
data <- rnorm(1000, mean = c(0, 4, 8), sd = c(1, 1.5, 0.5))

# Fit Gaussian mixture model
fit <- Mclust(data, G = 1:3)

# Extract optimal group number, means, and variances
optimal_G <- fit$G
optimal_means <- fit$parameters$mean
optimal_variances <- fit$parameters$variance

# Print the results
cat("Optimal Group Number:", optimal_G, "\n")
cat("Optimal Means:", "\n")
print(optimal_means)
cat("Optimal Variances:", "\n")
print(optimal_variances)

summary(fit)