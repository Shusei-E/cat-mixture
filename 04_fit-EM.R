library(tidyverse)
library(foreach)
library(glue)

data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")



# initialize mu {K x D x L}
init_mu_pr = matrix(rep(apply(data$y, 2, mean), data$K),
                     nrow = data$K,
                     byrow = TRUE)

mu <- init_mu <- array(NA, dim = c(data$K, data$D, data$L + 1))
# introduce some noise so clusters are not equivalent
noise = rnorm(data$K*data$D, sd = 0.05)

# y = 1 corresponds to second array
init_mu[, , data$L + 1] = (init_mu_pr - noise)
# y = 0 corresponds to first matrix
init_mu[, , data$L] = 1 - (init_mu_pr - noise)

# initialize theta {K x 1}
init_theta = rep(1/data$K, data$K)
init_theta = params$theta # correct thetas
theta = rep(NA, data$K)

# container for zeta {N x K}
zeta_hat = matrix(NA, nrow = data$N, ncol = data$K)


# iterations
iter = 1
store_iter = list()

while (iter < 100) {
  if (iter == 1) {
    theta = init_theta
    mu    = init_mu
  }

  # E step
  for (i in 1:data$N) {
    # responsibility for each k
    resp_i = rep(NA, data$K)

    for (k in 1:data$K) {
      # zeta_{ik} = theta_{k}prod^{D}_{j=1}prod^{L}_{l=0})(mu[k, j, (l + 1)]^(y[i, j] == l))
      resp_i[k] = foreach(j = 1:data$D, .combine = "*") %:%
        foreach(l = 0:(data$L - 1), .combine = "*") %do% {
          mu[k, j, (l + 1)]^(data$y[i, j] == l)
        }
    }

    numer_i = theta * resp_i # K x 1
    denom_i_k = sum(theta %*% resp_i) # scalar

    zeta_i = numer_i / denom_i_k # K x 1
    zeta_hat[i, ] = zeta_i
  }

  # M step

  for (k in 1:data$K) {
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

# check convergence

