library(rstan)
library(tidyverse)

data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")


# run sampling model
stanmodel <- stan_model("finite-mixture_stan-05-03.stan")
stanout <- rstan::sampling(stanmodel,
                           data,
                           chains = 1,
                           iter = 1000,
                           thin = 1,
													 )

stanmodel2 <- stan_model("finite-mixture_stan_bern.stan")
stanout <- rstan::vb(stanmodel2, data)

# save mixture model
write_rds(stanout, "data/mcmc/stan-fit.Rds")


# plot output --
theta <- params$theta
mu_true <- params$mu
pi_label <- glue("[{str_c(round(theta[1, ], 2), collapse = ', ')}]")
mu_label <- glue("[{str_c(round(mu_true[[1]], 2), collapse = ', ')}]")


stan_plot(stanout, "theta") +
  labs(caption = glue("True simplex is {pi_label}"))
ggsave("figures/sim_posterior_pi.png", w = 5, h = 3)

stan_trace(stanout, "theta", alpha = 0.3, nrow = 1)  +
  labs(caption = glue("True simplex is {pi_label}"))
ggsave("figures/sim_trace_pi.png", w = 10, h = 3)

stan_trace(stanout, glue("mu[1, {1:D}]"), alpha = 0.3, nrow = 2) +
  labs(caption = glue("True probabilities are {mu_label}"))
ggsave("figures/sim_trace_mu_cluster-01.png", w = 10, h = 3)



