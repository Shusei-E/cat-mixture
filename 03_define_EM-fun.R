library(tidyverse)
library(foreach)
library(glue)


# Diagnosis ----

# check sup-norm of changes in model parameters

#' Observed log likelihood
#'
#' @param t iteration index
#' @param obj the EM object with all iterations
#' @return the log likelihood
#'
loglik_obs <- function(t, obj = store_iter) {
  loglik_obs <- c()

  theta_t <- obj[[t]]$theta
  mu_t <- obj[[t]]$mu
  user_K <- length(theta_t)

  # for each N
  for (u in 1:data$U) {
    resp_i <- c()
    for (k in 1:user_K) {
      resp_i[k] = foreach(j = 1:data$D, .combine = "*") %:%
        foreach(l = 0:(data$L), .combine = "*") %do% {
          log(mu_t[k, j, (l + 1)])^(data$uy[u, j] == l)
        }
    }
    loglik_obs[u] <- data$n_u[u] * log(theta_t %*% resp_i)
  }
  sum(loglik_obs)
}



#'extract parameters from iteration t
#' @param t iteration index
#' @param obj the EM object with all iterations
#'
#' @return A tibble
vector_params <- function(t, loglik = FALSE, obj = store_iter) {
  theta_vector <- obj[[t]]$theta
  mu_vector <- obj[[t]]$mu[, , (data$L + 1)]
  df_i <- tibble(param_id = 1:(length(theta_vector) + length(mu_vector)),
                 type = c(rep("theta", length(theta_vector)),
                          rep("mu",    length(mu_vector))),
                 values = c(theta_vector, mu_vector)
  )
  if (loglik) df_i$loglik_obs <- loglik_obs(t)
  return(df_i)
}


#' Compute EM
cat_mixture <- function(data, user_K = 3, n_iter = 100, fast = TRUE) {

  # Unique profile and speed up? ----
  if (!fast) {
    data$U <- data$N
    data$n_u <- rep(1, data$N)
    data$uy <- data$y
  }

  # Setup -------
  # initialize theta {K x 1}
  init_theta = rep(1/user_K, user_K)
  # init_theta = params$theta # chat by giving it  correct thetas?

  # initialize mu {K x D x L}x
  init_Z_table = rmultinom(data$N, size = 1, prob = init_theta)
  init_Z = map_dbl(1:data$N, ~which(init_Z_table[, .x] == 1))

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


# Summarize ----

# stack all pre-post comparisons
summ_params <- function(store_iter, calc_loglik = TRUE) {
  foreach(t = 2:length(store_iter), .combine = "bind_rows") %do% {
    params_t <- vector_params(t, loglik = calc_loglik)
    params_tminus1 <- vector_params(t - 1, loglik = FALSE)

    left_join(params_tminus1, params_t,
              by = c("param_id", "type"),
              suffix = c("_pre", "_now")) %>%
      mutate(iter = t) %>%
      mutate(diff = abs(values_now - values_pre))
  }
}