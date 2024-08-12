# Load necessary library
library(MASS)  # For mvrnorm function to generate multivariate normal distributions

# Set the seed for reproducibility
set.seed(123)

# Define the dimensions
num_rows <- 500
num_cols <- 800

# Define the probabilities for the number of components
component_probs <- c(0.9, 0.09, 0.01)

# Function to generate a row from a Gaussian mixture
generate_gaussian_mixture_row <- function(num_cols) {
  # Randomly draw the number of components based on the given probabilities
  num_components <- sample(1:3, size = 1, prob = component_probs)
  
  # Initialize the row with zeros
  row <- numeric(num_cols)
  
  # Generate the components and their means and standard deviations
  for (i in 1:num_components) {
    # Determine the range for this component
    start_idx <- ifelse(i == 1, 1, floor((i - 1) * num_cols / num_components) + 1)
    end_idx <- floor(i * num_cols / num_components)
    
    # Generate mean and standard deviation for this component
    mu <- runif(1, min = 0.05, max = 0.95)
    sigma <- runif(1, min = 0.001, max = 0.05)
    
    # Populate the row for this component
    component_length <- end_idx - start_idx + 1
    row[start_idx:end_idx] <- mvrnorm(n = component_length, mu = mu, Sigma = matrix(sigma^2, nrow = 1))
  }
  
  # Ensure all values are between zero and one
  row <- pmax(0, pmin(1, row))
  return(row)
}

# Generate the matrix
gaussian_mixture_matrix <- t(replicate(num_rows, generate_gaussian_mixture_row(num_cols)))
