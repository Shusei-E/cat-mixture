library(purrr)
library(rstan)
library(brms)
library(glue)

set.seed(02138)

# dimensions
L <- 1L
D <- 8L
N <- 300L
K <- 3L

# hyperparameter
alpha <- (K:1)


# cluster assignment
theta <- rdirichlet(1, alpha)
Z_table <- rmultinom(N, 1, theta)
Z <-  map_dbl(1:N,  ~which(Z_table[, .x] == 1))

# set theta parameters
mu <- list(
  `1` = rep(0.05, D),
  `2` = rep(0.95, D),
  `3` = rbeta(D, 2, 5)
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
as_tibble(k_vanilla$centers) %>%
  mutate(n = k_vanilla$size)

# check sample mean | correct cluster assignment
as_tibble(data$y) %>%
  mutate(cluster = Z) %>%
  group_by(cluster) %>%
  summarize_all(mean)

as_tibble(data$y) %>%
  mutate(profile = glue("{V1}{V2}{V3}{V4}{V5}{V6}{V7}{V8}")) %>%
  count(profile, sort = TRUE)
