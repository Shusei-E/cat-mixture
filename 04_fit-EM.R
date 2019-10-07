library(tidyverse)
library(foreach)
library(glue)

# Data --------
data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")


# Setup -------
user_K <- 2

# initialize theta {K x 1}
init_theta = rep(1/user_K, user_K)
# init_theta = params$theta # chat by giving it  correct thetas?

# initialize mu {K x D x L} -- currently only supports L = 2
init_Z_table <- rmultinom(data$N, size = 1, prob = init_theta)
init_Z <- map_dbl(1:data$N, ~which(init_Z_table[, .x] == 1))
init_muhat <- flatten_dbl(map(1:user_K, ~colMeans(data$y[init_Z == .x, ])))

init_mu_pr = matrix(init_muhat,
                    nrow = user_K,
                    byrow = TRUE)

mu <- init_mu <- array(NA, dim = c(user_K, data$D, data$L + 1))

# y = 1 corresponds to second array
init_mu[, , data$L + 1] = init_mu_pr
# y = 0 corresponds to first matrix
init_mu[, , data$L] = 1 - init_mu_pr

# container for theta
theta = rep(NA, user_K)

# container for zeta {N x K}
zeta_hat = matrix(NA, nrow = data$N, ncol = user_K)





# fit EM -----------
# iterations
iter = 1
store_iter = list()

while (iter <= 50) {
  if (iter == 1) {
    theta = init_theta
    mu    = init_mu
  }

  # E step
  for (i in 1:data$N) {
    # responsibility for each k
    resp_i = rep(NA, user_K)

    for (k in 1:user_K) {
      # zeta_{ik} = theta_{k}prod^{D}_{j=1}prod^{L}_{l=0})(mu[k, j, (l + 1)]^(y[i, j] == l))
      resp_i[k] = foreach(j = 1:data$D, .combine = "*") %:%
        foreach(l = 0:(data$L), .combine = "*") %do% {
          mu[k, j, (l + 1)]^(data$y[i, j] == l)
        }
    }

    numer_i = theta * resp_i # K x 1
    denom_i_k = sum(theta %*% resp_i) # scalar

    zeta_i = numer_i / denom_i_k # K x 1
    zeta_hat[i, ] = zeta_i
  }

  # M step

  for (k in 1:user_K) {
    sum_zeta_k = sum(zeta_hat[, k])

    # update theta
    theta[k] = (1/data$N)*sum_zeta_k

    # update mu
    for (l in 0:(data$L - 1)) {
      y_matches_l = data$y[, j] == l
      mu[k, j, (l + 1)]  = sum(y_matches_l %*% zeta_hat[, k]) / sum_zeta_k
    }
    # last category is 1 minus the rest
    mu[k, j, (data$L + 1)] =  1 - sum(mu[k, j, 1:(data$L)])
  }

  # store each iter
  store_iter[[iter]] = list(mu = mu, theta = theta, zeta = zeta_hat)
  cat(glue("iter: {iter}"), "\n")

  iter <- iter + 1
}



# Store -------
write_rds(store_iter, "data/EM/sim-iterations.Rds")
