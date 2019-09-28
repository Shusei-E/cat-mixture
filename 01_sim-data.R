library(purrr)
library(rstan)
library(brms)

# dimensions
M <- 3L
K <- 5L
J <- 10L
N <- 400L

# hyperparameter
alpha <- c(2.0, 1.5, 1.0)

# priors
psi <- rdirichlet(1, alpha)
Z_table <- rmultinom(N, 1, psi)
Z <-  map_dbl(1:N,  ~which(Z_table[, .x] == 1) - 1)

# set parameters
theta_z <- list(
  `0` = c(0.97, 0.01, 0.02),
  `1` = c(0.01, 0.97, 0.02),
  `2` = c(0.49, 0.49, 0.02)
)

# Generate data
Y <- array(NA, dim = c(N, J))
for (j in 1:J) {
  theta_i <- map(as.character(Z), ~ theta_z[[.x]])
  y_j <- map_dbl(theta_i, ~which(rmultinom(1, 1, .x) == 1) - 1)
  Y[, j] <- y_j
}

