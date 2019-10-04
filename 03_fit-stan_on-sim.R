library(rstan)
library(tidyverse)

data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")


# run sampling model
stanmodel <- stan_model("finite-mixture_stan-05-03.stan")
stanout <- rstan::sampling(stanmodel,
                           data,
                           chains = 4,
                           iter = 2e3,
                           thin = 1)

# save mixture model
write_rds(stanout, "data/mcmc/stan-fit.Rds")


# plot output --
pi_mixture <- params$pi
mu_true <- params$mu
pi_label <- glue("[{str_c(round(pi_mixture[1, ], 2), collapse = ', ')}]")
mu_label <- glue("[{str_c(round(mu_true[[1]], 2), collapse = ', ')}]")


stan_plot(stanout, "pi_mixture") +
  labs(caption = glue("True simplex is {pi_label}"))
ggsave("figures/sim_posterior_pi.png", w = 5, h = 3)

stan_trace(stanout, "pi_mixture", alpha = 0.3, nrow = 1)  +
  labs(caption = glue("True simplex is {pi_label}"))
ggsave("figures/sim_trace_pi.png", w = 10, h = 3)

stan_trace(stanout, glue("mu[1, {1:D}]"), alpha = 0.3, nrow = 2) +
  labs(caption = glue("True probabilities are {mu_label}"))
ggsave("figures/sim_trace_mu_cluster-01.png", w = 10, h = 3)



