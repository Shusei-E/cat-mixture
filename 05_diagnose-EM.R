library(tidyverse)
library(foreach)
library(ggthemes)
library(lemon)
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
  for (u in 1:data$U) {
    resp_i <- c()
    for (k in 1:user_K) {
      resp_i[k] = foreach(j = 1:data$D, .combine = "+") %:%
        foreach(l = 0:(data$L), .combine = "+") %do% {
          (data$uy[u, j] == l)*log(mu_t[k, j, (l + 1)])
        }
    }
    loglik_obs[u] <- data$n_u[u] * prod(log(theta_t) + resp_i[1:user_K])
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

# stack all pre-post comparisons
params_stacked <- foreach(t = 2:length(store_iter), .combine = "bind_rows") %do% {
  params_t <- vector_params(t, loglik = TRUE)
  params_tminus1 <- vector_params(t - 1, loglik = FALSE)

  left_join(params_tminus1, params_t,
            by = c("param_id", "type"),
            suffix = c("_pre", "_now")) %>%
    mutate(iter = t) %>%
    mutate(diff = abs(values_now - values_pre))
}


# plot max of diff
params_stacked %>%
  mutate(llobs_scale = loglik_obs / (data$N*data$D)) %>%
  group_by(iter) %>%
  summarize(`Maximum Change in Parameter (probability scale)` = max(diff),
            `Observed Log Likelihood (per data point)` = unique(llobs_scale)) %>%
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
ggsave("figures/sim_EM_change-in-params.pdf", w = 4, h = 4)
