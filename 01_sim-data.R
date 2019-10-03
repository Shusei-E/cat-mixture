library(purrr)
library(rstan)
library(brms)
library(glue)

set.seed(02138)

# dimensions
M <- 2L
K <- 5L
D <- 8L
N <- 400L

# hyperparameter
alpha <- c(3, 1, 1, 1, 1)


# cluster assignment
pi <- rdirichlet(1, alpha)
Z_table <- rmultinom(N, 1, pi)
Z <-  map_dbl(1:N,  ~which(Z_table[, .x] == 1))

# set parameters
mu <- list(
  `1` = rep(.01, D),
  `2` = rep(.05, D),
  `3` = rep(.50, D),
  `4` = rep(.99, D),
  `5` = rbeta(D, 1, 1)
)

# each mu vector is K by D
stopifnot(all(map_lgl(mu, ~length(.x) == D)))


# Generate data
y <- array(NA, dim = c(N, D))
for (i in 1:N) {
  y[i, ] <- rbinom(D, size = 1, prob = mu[[Z[i]]])
}

# put together data
data <- list(D = D,
             K = K,
             N = N,
             y = y)

# target params
params <- list(pi = pi,
               mu = mu,
               Z = Z)

write_rds(data, "data/sim-data.Rds")
write_rds(params, "data/sim-params.Rds")

# check vanilla k means
k_vanilla <- kmeans(data$y, centers = K)
kdf_vanilla <- as_tibble(k_vanilla$centers) %>%
  mutate(n = k_vanilla$size)