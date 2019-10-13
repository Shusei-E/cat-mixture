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


# cluster assignment
theta <- rdirichlet(1, alpha)
Z_table <- rmultinom(N, 1, theta)
Z <-  map_dbl(1:N,  ~which(Z_table[, .x] == 1))

# set theta parameters
mu <- list(
  `1` = rep(0.10, D),
  `2` = rep(0.90, D),
  `3` = rbeta(D, 2, 5)
)

# each mu vector is K by D
stopifnot(all(map_lgl(mu, ~length(.x) == D)))


# Generate data
y <- array(NA, dim = c(N, D))
for (i in 1:N) {
  y[i, ] <- rbinom(D, size = L, prob = mu[[Z[i]]])
}

# unique profiles
unique_y <- as_tibble(data$y) %>%
  mutate(voter = 1:n()) %>%
  pivot_longer(-voter, names_to = "j", values_to = "y") %>%
  group_by(voter) %>%
  summarize(profile = str_c(y, collapse = "")) %>%
  count(profile, name = "n_u") %>%
  separate(profile, into = str_c("v", 1:D), sep = 1:D) %>%
  mutate_all(as.integer)

unique_y_mat <- unique_y[, 1:D]
n_u <- unique_y$n_u
U <- length(n_u)


# put together data
data <- list(D = D,
             K = K,
             N = N,
             L = L,
             y = y,
             uy = unique_y_mat,
             n_u = n_u,
             U = U,
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
