library(tidyverse)
library(foreach)
library(glue)

#' Compute EM
cat_mixture <- function(data, user_K = 3, n_iter = 100, fast = TRUE) {

  # Unique profile and speed up? ----
  if (!fast) {
    data$U <- data$N
    data$n_u <- rep(1, data$N)
    data$uy <- data$y
  }

  # Setup -------
  user_K <- 3

  # initialize theta {K x 1}
  init_theta = rep(1/user_K, user_K)
  # init_theta = params$theta # chat by giving it  correct thetas?

  # initialize mu {K x D x L}x
  init_Z_table <- rmultinom(data$N, size = 1, prob = init_theta)
  init_Z <- map_dbl(1:data$N, ~which(init_Z_table[, .x] == 1))

  # data$L
  mu <- init_mu <- array(NA, dim = c(user_K, data$D, data$L + 1))
  for (l in 0:data$L) {
    mu_vector <- flatten_dbl(map(1:user_K, ~colMeans(data$y[init_Z == .x, ] == l)))
    init_mu[, , l + 1] = matrix(mu_vector, nrow = user_K, byrow = TRUE)
  }

  # container for theta
  theta = rep(NA, user_K)

  # container for zeta {N x K}
  zeta_hat = matrix(NA, nrow = data$U, ncol = user_K)



  # iterations ------
  iter = 1
  store_iter = list()

  # start EM loop ------
  while (iter <= n_iter) {
    if (iter == 1) {
      theta = init_theta
      mu    = init_mu
    }

    # E step
    for (u in 1:data$U) {
      # responsibility of type k
      resp_u = rep(NA, user_K)

      for (k in 1:user_K) {
        # zeta_{ik} = theta_{k}prod^{D}_{j=1}prod^{L}_{l=0})(mu[k, j, (l + 1)]^(y[i, j] == l))
        resp_u[k] = foreach(j = 1:data$D, .combine = "*") %:%
          foreach(l = 0:(data$L), .combine = "*") %do% {
            (mu[k, j, (l + 1)])^(data$uy[u, j] == l)
          }
      }

      numer_u = theta * resp_u # K x 1
      denom_u_k = sum(theta %*% resp_u) # scalar

      zeta_u = numer_u / denom_u_k # K x 1
      zeta_hat[u, ] = zeta_u
    }

    # M step

    for (k in 1:user_K) {
      sum_zeta_k = data$n_u %*% zeta_hat[, k]

      # update theta
      theta[k] = (1/data$N)*sum_zeta_k
      # correct theta
      # theta[k] = params$theta[k]

      # update mu
      for (j in 1:data$D) {
        for (l in data$L:1) {
          # 1(Y_{ij} = ell)
          y_matches_l = (data$uy[, j] == l)
          # update mu
          mu[k, j, (l + 1)]  = sum(data$n_u * y_matches_l * zeta_hat[, k]) /
            sum_zeta_k
        }

        # last category is 1 minus the rest
        mu[k, j, (0 + 1)] =  1 - sum(mu[k, j, ((1:data$L) + 1)])
      }

    }

    # store each iter
    store_iter[[iter]] = list(mu = mu, theta = theta, zeta = zeta_hat)
    cat(glue("iter: {iter}"), "\n")

    iter = iter + 1
  }
  return(store_iter)
}