# Cpp
Rcpp::sourceCpp("02_cpp-functions.cpp")

# Diagnosis ----

# check sup-norm of changes in model parameters

#' Observed log likelihood
#'
#' @param t iteration index
#' @param obj the EM object with all iterations
#' @param data the dataset with observed values
#' @param fast logical, if the object was computed
#'  by collpasing to profiles (TRUE) or not (FALSE)
#'
#'
#' @return the log likelihood
#'
loglik_obs <- function(t, obj, data, fast = TRUE, IIA = FALSE) {
  loglik_obs <- c()

  if (IIA) fast <- FALSE

  if (!fast) {
    data$U <- data$N
    data$n_u <- rep(1, data$N)
    data$uy <- data$y
  }

  theta_t <- obj[[t]]$theta
  mu_t <- obj[[t]]$mu
  user_K <- length(theta_t)

  # for each N
  for (u in 1:data$U) {
    resp_i = c()
    for (k in 1:user_K) {
      # if not IIA, then follow original
      if (!IIA) {
        resp_i[k] = mu_yvec(mu_k = mu_t[k, , ], y = data$uy[u, ], L = data$L)
      }
      # recale if IIA
      if (IIA) {
        resp_i_k = c()
        for (j in 1:data$D) {
          S = switch(data$m[u, j], `1` = c(0, 1), `2` = c(0, 2), `3` = c(0, 1, 2))
          sum_mu = sum(mu_t[k, j, S + 1])
          resp_i_k[j] = foreach(l = S, .combine = "*") %do% {
            (mu_t[k, j, (l + 1)] / sum_mu)^(data$uy[u, j] == l)
          }
        }
        resp_i[k] = prod(resp_i_k)
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
vector_params <- function(t, loglik = FALSE, obj, data, IIA = FALSE) {
  theta_vector <- obj[[t]]$theta
  mu_vector    <- obj[[t]]$mu[, , (1:(data$L + 1))]

  user_K <- length(theta_vector)
  theta_names <- str_c("theta_", 1:user_K)
  mu_names    <- str_c(str_c(str_c("mu_", 1:user_K, "_"),
                             rep(1:data$D, each = user_K), "_"),
                       rep(0:data$L, each = user_K*data$D))

  df_i <- tibble(param_id = c(theta_names, mu_names),
                 values = c(theta_vector, mu_vector)
  )
  if (loglik) df_i$loglik_obs <- loglik_obs(t, obj, data, IIA = IIA)
  return(df_i)
}


#' Compute EM
cat_mixture <- function(data, user_K = 3, n_iter = 100, fast = TRUE, IIA = FALSE,
                        init = c("kmeans", "equal"),
                        seed = 02138,
                        theta = NULL, mu = NULL, zeta_hat  = NULL) {

  # checkpackages
  stopifnot(all(require(purrr), require(foreach), require(glue)))
  if (IIA) stopifnot(require(mlogit))

  # cannot handle groupings in IIA case
  if (IIA) fast <- FALSE

  # Unique profile and speed up? ----
  if (!fast) {
    data$U <- data$N
    data$n_u <- rep(1, data$N)
    data$uy <- data$y
  }

  # Setup and initialize -------
  # unless all values are given, start from initial guesses
  if (any(is.null(theta), is.null(mu), is.null(zeta_hat))) {

    if (init == "kmeans") {
    # run k-means on binarized data
    init_binary <- matrix(recode(data$y, `2` = 1, `1` = 0,`0` = 0),
                          nrow = data$N,
                          ncol = data$D,
                          byrow = FALSE)

    set.seed(seed)
    k_init <- kmeans(init_binary, centers = user_K)
    # pre-sort so largest cluster tends to be at front
    c_order <- order(table(k_init$cluster), decreasing = TRUE)

    switcher <- 1:user_K
    names(switcher) <- c_order

    k_cl <- dplyr::recode(k_init$cluster, !!!switcher)

    # initialize theta {K x 1}
    init_theta = prop.table(table(k_init$cluster))[c_order]

    # initialize mu {K x D x L}x
    init_Z = k_cl

    # data$L
    mu <- init_mu <- array(NA, dim = c(user_K, data$D, data$L + 1))
    for (l in 0:data$L) {
      mu_vector <- flatten_dbl(map(1:user_K, ~colMeans(data$y[init_Z == .x, ] == l)))
      init_mu[, , l + 1] = matrix(mu_vector, nrow = user_K, byrow = TRUE)
    }

    # correct so that no cell is exactly zero
    for (j in 1:data$D) {
      for (k in 1:user_K) {
        if (any(init_mu[k, j, ] < 1e-5))
          init_mu[k, j, ] = (init_mu[k, j, ]  + 0.05) /  sum((init_mu[k, j, ]  + 0.05))
      }
    }

    # container for theta
    theta = rep(NA, user_K)

    # container for zeta {N x K}
    zeta_hat = matrix(NA, nrow = data$U, ncol = user_K)

    # initialize
    theta = init_theta
    mu    = init_mu
    }

    if (init == "equal") {
      # initialize theta {K x 1}
      init_theta = rep(1/user_K, user_K)

      # initialize mu {K x D x L}
      set.seed(seed)
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

      # initialize
      theta = init_theta
      mu    = init_mu
    }
  }

  # iterations ------
  iter = 1
  store_iter = list()

  # start EM loop ------
  while (iter <= n_iter) {
    ## E step ----
    t_e_start <- Sys.time()
    for (u in 1:data$U) {
      # responsibility of type k
      resp_u = rep(NA, user_K)

      for (k in 1:user_K) {
        if (!IIA) {
          resp_u[k] = mu_yvec(mu_k = mu[k, , ], y = data$uy[u, ], L = data$L)
        }

        if (IIA) {
          stopifnot(data$L == 2)
          resp_u_k <- c()
          for (j in 1:data$D) {
            m = data$m[u, j]
            S_m = switch(m, `1` = c(0, 1), `2` = c(0, 2), `3` = c(0, 1, 2))
            sum_mu = sum(mu[k, j, (S_m + 1)])

            resp_u_k[j] <- foreach(l = S_m, .combine = "*") %do% {
              (mu[k, j, (l + 1)] / sum_mu)^(data$uy[u, j] == l)
            }
          }
          resp_u[k] = prod(resp_u_k)
        }
      }

      numer_u = theta * resp_u # K x 1
      denom_u_k = sum(theta %*% resp_u) # scalar

      zeta_u = numer_u / denom_u_k # K x 1
      zeta_hat[u, ] = zeta_u
    }
    t_e_end <- Sys.time()

    # M step ----
    t_m_start <- Sys.time()
    for (k in 1:user_K) {
      sum_zeta_k = data$n_u %*% zeta_hat[, k]

      # update theta
      theta[k] = (1/data$N)*sum_zeta_k
      # correct theta
      # theta[k] = params$theta[k]

      # update mu
      for (j in 1:data$D) {

        # default
        if (!IIA) {
          for (l in data$L:1) {
            # 1(Y_{ij} = ell)
            y_matches_l = (data$uy[, j] == l)

            # update mu
            mu[k, j, (l + 1)]  = sum(data$n_u * y_matches_l * zeta_hat[, k]) /
              sum_zeta_k
          }
          # last category is 1 minus the rest
          mu[k, j, (0 + 1)] =  1 - sum(mu[k, j, ((data$L:1) + 1)])
        }

        # IIA model
        if (IIA) {
          # regression data with missingnesses dropped
          y_no_na <- tibble(id   = rep(1:data$U, each = data$L + 1),
                            y_j  = rep(data$uy[, j], each = data$L + 1),
                            M_j  = rep(data$m[, j], each = data$L + 1),
                            zeta = rep(zeta_hat[, k], each = (data$L + 1)),
                            alt  = rep(0:data$L, data$U)) %>%
            mutate(y_chosen = y_j == alt) %>%
            filter(!(M_j == 1 & alt == 2), !(M_j == 2 & alt == 1)) %>%
            as.data.frame()

          mlogit_d <- mlogit.data(y_no_na,
                                  choice  = "y_chosen",
                                  chid.var = "id",
                                  alt.var = "alt",
                                  shape = "long")

          # multinomial logit
          mfit <- mlogit(y_chosen ~ 1, weights = zeta, mlogit_d)

          # rescale to probabilities
          mfit_coefs_e <- exp(c(0, coef(mfit)))
          mu[k, j, ] <- mfit_coefs_e / sum(mfit_coefs_e)
        }
      }
    }
    t_m_end <- Sys.time()

    # store each iter
    opts <- list(fast = fast,
                 IIA = IIA,
                 N = data$N,
                 D = data$D,
                 init = init,
                 time_e = difftime(t_e_end, t_e_start),
                 time_m = difftime(t_m_end, t_m_start))
    store_iter[[iter]] = list(mu = mu, theta = theta, zeta = zeta_hat, opts = opts)
    cat(glue("iter: {iter}"), "\n")

    iter = iter + 1
  }
  return(store_iter)
}


# Summarize ----

# stack all pre-post comparisons
summ_params <- function(store_iter, data, calc_loglik = TRUE, IIA = FALSE) {
  foreach(t = 2:length(store_iter), .combine = "bind_rows") %do% {
    t_loglik_start <- Sys.time()
    params_t       <- vector_params(t, loglik = calc_loglik, store_iter, data, IIA = IIA)
    t_loglik_end <- Sys.time()
    params_tminus1 <- vector_params(t - 1, loglik = FALSE, store_iter, data, IIA = IIA)

    left_join(params_tminus1, params_t,
              by = c("param_id"),
              suffix = c("_pre", "_now")) %>%
      mutate(iter = t,
             diff = abs(values_now - values_pre),
             time_l = difftime(t_loglik_end, t_loglik_start))
  }
}


#' graph parameter fit
trend_stacked <- function(params_stacked, dat = data) {

  if (!"loglik_obs" %in% colnames(params_stacked)) {
    summ_df <-  group_by(params_stacked, iter) %>%
      summarize(`Maximum Change in Parameter (probability scale)` = max(diff))
  }

  if ("loglik_obs" %in% colnames(params_stacked)) {
    summ_df <- params_stacked %>%
      mutate(llobs_scale = loglik_obs / (dat$N*dat$D)) %>%
      group_by(iter) %>%
      summarize(`Maximum Change in Parameter (probability scale)` = max(diff),
                `Observed Log Likelihood (per data point)` = unique(llobs_scale))
  }

  summ_df %>%
    pivot_longer(cols = -c(iter), names_to = "metric", values_to = "value") %>%
    ggplot(aes(iter, value)) +
    facet_rep_wrap(~metric, scales = "free_y", ncol = 1) +
    coord_capped_cart(bottom='both', left = 'both') +
    geom_point(size = 0.5) +
    geom_line() +
    theme_clean() +
    theme(plot.background = element_rect(color = NA),
          axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.caption = element_text(size = 6),
          strip.background = element_rect(fill = "lightgray")) +
    labs(x = "EM Iteration",
         y = "Metric",
         caption = glue("Note: The parmater vector is the estimated theta's and mu's combined.
                      Higher observed log likelihood and lower sup-norms both indicate better fit."))
}

