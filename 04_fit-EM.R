library(tidyverse)
library(foreach)
library(glue)

data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")



# initialize mu {K x D x L}
init_mu_pr <- matrix(rep(apply(data$y, 2, mean), data$K),
                     nrow = data$K,
                     byrow = TRUE)
mu <- init_mu <- array(NA, dim = c(data$K, data$D, data$L + 1))
init_mu[, , data$L + 1] <- init_mu_pr # y = 1 corresponds to second array
init_mu[, , data$L] <- 1 - init_mu_pr # y = 0 corresponds to first matrix

# initialize theta {K x 1}
init_theta <- rep(1/data$K, data$K)
theta <- rep(NA, data$K)

# container for zeta
zeta_hat <- matrix(NA, nrow = data$N, ncol = data$K)

iter <- 1

while (iter < 10) {
  if (iter == 1) {
    theta <- init_theta
    mu <- init_mu
  }

  # E step
  for (i in 1:data$N) {
    # responsibility for each k
    resp_i_k <- rep(NA, data$K)


    for (k in 1:data$K) {
      # zeta_{ik} = theta_{k}prod^{D}_{j=1}prod^{L}_{l=0})(mu[k, j, (l + 1)]^(y[i, j] == l))
      resp_i_k[k] <- foreach(j = 1:data$D, .combine = "*") %:%
        foreach(l = 0:(data$L - 1), .combine = "*") %do% {
          mu[k, j, (l + 1)]^(data$y[i, j] == l)
        }
    }

    numer_i <- theta * resp_i_k # K x 1
    denom_i_k <- sum(theta %*% resp_i_k) # scalar

    zeta_i <- numer_i / denom_i_k # K x 1
    zeta_hat[i, ] <- zeta_i
  }

  # M step

  for (k in 1:data$K) {
    sum_zeta_k <- sum(zeta_hat[, k])

    theta[k] = (1/data$N)*sum_zeta_k

    for (l in 0:data$L) {
      y_matches_l = data$y[, j] == l
      mu[k, j, (l + 1)]  = sum(y_matches_l %*% zeta_hat[, k]) / sum_zeta_k
    }
  }
  iter <- iter + 1
  cat(glue("iter: {iter}"), "\n")
}

# check convergence

