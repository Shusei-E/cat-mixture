library(tidyverse)
library(foreach)
library(glue)

# Data --------
data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")

store_iter <- read_rds("data/EM/sim-iterations.Rds")

# Diagnostics ---------

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
  for (i in 1:data$N) {
    resp_i <- c()
    for (k in 1:user_K) {
      resp_i[k] = foreach(j = 1:data$D, .combine = "+") %:%
        foreach(l = 0:(data$L), .combine = "+") %do% {
          (data$y[i, j] == l)*log(mu_t[k, j, (l + 1)])
        }
    }
    loglik_obs[i] <- prod(log(theta_t) + resp_i[1:user_K])
  }
  sum(loglik_obs)
}

loglik_obs(1)
loglik_obs(10)

#'extract parameters from iteration t
#' @param t iteration index
#' @param obj the EM object with all iterations
#'
#' @return A tibble
vector_params <- function(t, obj = store_iter) {
  theta_vector <- obj[[t]]$theta
  mu_vector <- obj[[t]]$mu[, , (data$L + 1)]
  df_i <- tibble(param_id = 1:(length(theta_vector) + length(mu_vector)),
                 type = c(rep("theta", length(theta_vector)),
                          rep("mu",    length(mu_vector))),
                 values = c(theta_vector, mu_vector)
  )
  return(df_i)
}

# stack all pre-post comparisons
params_stacked <- foreach(t = 2:length(store_iter), .combine = "bind_rows") %do% {
  params_t <- vector_params(t)
  params_tminus1 <- vector_params(t - 1)

  left_join(params_tminus1, params_t,
            by = c("param_id", "type"),
            suffix = c("_pre", "_now")) %>%
    mutate(iter = t) %>%
    mutate(diff = values_now - values_pre)
}


# plot max of diff
params_stacked %>%
  group_by(iter) %>%
  summarize(max_diff = max(diff)) %>%
  ggplot(aes(iter, max_diff)) +
  geom_point() +
  geom_line() +
  labs(x = "Iteration",
       y = "max(change in parameter value)",
       caption = glue("the parmater vector is {data$K} theta's and {data$K*data$D} mu's combined."))
ggsave("figures/sim_EM_change-in-params.png", w = 4, h = 3)
