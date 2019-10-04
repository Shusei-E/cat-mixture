library(purrr)
library(rstan)
library(brms)
library(glue)

set.seed(02138)

# dimensions
M <- 2L
K <- 3L
D <- 5L
N <- 400L

# hyperparameter
alpha <- c(1, 1, 1)


# cluster assignment
theta <- rdirichlet(1, alpha)
Z_table <- rmultinom(N, 1, theta)
Z <-  map_dbl(1:N,  ~which(Z_table[, .x] == 1))

# set parameters
mu <- list(
  `1` = rnorm(D, 0, 3),
  `2` = rnorm(D, 0, 3),
  `3` = rnorm(D, 0, 3)
  # `4` = rnorm(D, 0, 3),
  # `5` = rnorm(D, 0, 3)
)

# each mu vector is K by D
stopifnot(all(map_lgl(mu, ~length(.x) == D)))


# Generate data
y <- array(NA, dim = c(N, D))
for (i in 1:N) {
  # y[i, ] <- rbinom(D, size = 1, prob = mu[[Z[i]]])
  y[i, ] <- rnorm(D, mean = mu[[Z[i]]], sd=0.05)
}

# put together data
data <- list(D = D,
             K = K,
             N = N,
             y = y,
             alpha = alpha)

# target params
params <- list(theta = theta,
               mu = mu,
               Z = Z)

write_rds(data, "data/sim-data.Rds")
write_rds(params, "data/sim-params.Rds")

# check vanilla k means
k_vanilla <- kmeans(data$y, centers = K)
kdf_vanilla <- as_tibble(k_vanilla$centers) %>%
  mutate(n = k_vanilla$size)

# Gaussian Mixture
library(mclust)
gmm <- Mclust(y, G=K)
summary(gmm)
table(Z)
