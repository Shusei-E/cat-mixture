library(rstan)

mod <- stan_model("sims/mc-stan_05-03.stan")

fit <- sampling(mod,
                data = list(K = 5, N = 100, y = runif(100)))
