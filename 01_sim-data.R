library(purrr)
library(rstan)
library(brms)
library(glue)

set.seed(02138)

# dimensions
L <- 2L
D <- 5L
N <- 300L
K <- 3L

# hyperparameter
alpha <- (K:1)^(1.5)

#  assume missing by not available data?
missing <- FALSE

# cluster assignment
theta <- rdirichlet(1, alpha)
Z_table <- rmultinom(N, 1, theta)
Z <-  map_dbl(1:N,  ~which(Z_table[, .x] == 1))

# set theta parameters
mu_k <- list(
  `1` = c(0.05, 0.05, 0.90),
  `2` = c(0.10, 0.10, 0.80),
  `3` = c(0.15, 0.20, 0.65)
)

table(Z)

# data$L
mu <- array(NA, dim = c(K, D, L + 1))
for (k in 1:K) {
  for (j in 1:D) {
    # scales vary by k primarily, but nosie by j
    mu[k, j, ] <- mu_k[[k]] + 0.01*k*(j - 1)
    # rescale
    mu[k, j, ] <- mu[k, j, ] / sum(mu[k, j, ])
  }
}

# Missingness
if (missing)
  m <- array(sample(1:3, size = N*D, replace = TRUE), dim = c(N, D))
if (!missing)
  m <- array(3, dim = c(N, D))

# Generate data
y <- array(NA, dim = c(N, D))
indiv_mu <- array(NA, dim = c(N, D, L + 1))
for (i in 1:N) {
  for (j in 1:D) {
    mu_i <- mu[Z[i], j, ]

    l_not_available <- switch(m[i, j], `1` = 2, `2` = 1, `3` = NA)
    if (!is.na(l_not_available)) mu_i[1 + l_not_available] <- 0

    y[i, j] <- sample(0:L, size = 1, prob = mu_i)
    indiv_mu[i, j, ] <- mu_i
  }
}


# unique profiles
unique_y <- as_tibble(data$y) %>%
  mutate(voter = 1:n()) %>%
  pivot_longer(-voter, names_to = "j", values_to = "y") %>%
  group_by(voter) %>%
  summarize(profile = str_c(y, collapse = "")) %>%
  count(profile, name = "n_u") %>%
  separate(profile, into = str_c("v", 1:D), sep = 1:D) %>%
  mutate_all(as.integer) %>%
  as.matrix()

unique_y_mat <- unique_y[, 1:D]
n_u <- unique_y[, "n_u"]
U <- length(n_u)


# put together data
data <- list(D = D,
             K = K,
             N = N,
             L = L,
             y = y,
             m = m,
             uy = unique_y_mat,
             n_u = n_u,
             U = U,
             alpha = alpha)

# target params
params <- list(theta = theta,
               mu = mu,
               Z = Z)

if (missing)
  write_rds(data, "data/sim-data_missing.Rds")

if (!missing)
  write_rds(data, "data/sim-data.Rds")

write_rds(params, "data/sim-params.Rds")

cbind(Z = params$Z, data$y) %>%
  as_tibble() %>%
  group_by(Z) %>%
  summarize_all(~mean(.x == 2))

table(params$Z)
