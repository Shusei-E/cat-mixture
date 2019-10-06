library(purrr)
library(rstan)
library(brms)
library(glue)

set.seed(021382)

# dimensions
M <- 2L
K <- 5L
D <- 8L
L <- 1L
N <- 100L

# hyperparameter
alpha <- (K:1)^2


# cluster assignment
theta <- rdirichlet(1, alpha)
Z_table <- rmultinom(N, 1, theta)
Z <-  map_dbl(1:N,  ~which(Z_table[, .x] == 1))

# setpi parameters
mu <- list(
  `1` = rep(.01, D),
  `2` = rep(.05, D),
  `3` = rep(.10, D),
  `4` = rep(.99, D),
  `5` = rbeta(D, 2, 5)
)

# each mu vector is K by D
stopifnot(all(map_lgl(mu, ~length(.x) == D)))


# Generate data
y <- array(NA, dim = c(N, D))
for (i in 1:N) {
  y[i, ] <- rbinom(D, size = L, prob = mu[[Z[i]]])
}

# put together data
data <- list(D = D,
             K = K,
             N = N,
             L = L,
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