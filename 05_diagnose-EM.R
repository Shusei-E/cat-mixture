library(tidyverse)
library(foreach)
library(glue)

# Data --------
data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")

store_iter <- read_rds("data/EM/sim-iterations.Rds")

# Diagnostics ---------

# check sup-norm of changes in model parameters

# extract parameters
vector_params <- function(i, obj = store_iter) {
  theta_vector <- obj[[i]]$theta
  mu_vector <- obj[[i]]$mu[, , (data$L + 1)]
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
