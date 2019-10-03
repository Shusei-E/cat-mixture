library(rstan)
library(tidyverse)

data <- read_rds("data/sim-data.Rds")

# run sampling model
stanmodel <- stan_model("finite-mixture_stan-05-03.stan")
stanout <- rstan::sampling(stanmodel,
                           data,
                           chains = 5,
                           iter = 1e3,
                           thin = 5)

# save mixture model
write_rds(stanout, "data/mcmc/stan-fit.Rds")

# plot output --
stan_plot(stanout, "pi")
ggsave("figures/sim_posterior_pi.png", w = 5, h = 3)

stan_trace(stanout, "pi", alpha = 0.5, nrow = 1)
ggsave("figures/sim_trace_pi.png", w = 8, h = 4)

stan_trace(stanout, glue("mu[1, {1:D}]"), alpha = 0.5, nrow = 2)
ggsave("figures/sim_trace_mu_cluster-01.png", w = 8, h = 4)



