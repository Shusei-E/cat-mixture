library(tidyverse)
library(foreach)
library(patchwork)
library(ggthemes)
library(lemon)
library(glue)

source("03_define_EM-fun.R")

# Data --------
data <- read_rds("data/sim-data.Rds")
data_miss <- read_rds("data/sim-data_missing.Rds")
params <- read_rds("data/sim-params.Rds")
store_iter_1 <- read_rds("data/EM/sim-iterations.Rds")
store_iter_2 <- read_rds("data/EM/sim-iterations_IIA-full.Rds")
store_iter_3 <- read_rds("data/EM/sim-iterations_IIA-miss.Rds")


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



# Diagnostics ---------
pst <- summ_params(store_iter, data = data, calc_loglik = TRUE)
pst_1 <- summ_params(store_iter_1, data = data, calc_loglik = TRUE)
pst_2 <- summ_params(store_iter_2, data = data, calc_loglik = TRUE)
pst_3 <- summ_params(store_iter_3, data = data_miss, calc_loglik = TRUE)


trend_stacked(pst)
trend_stacked(pst_1) + trend_stacked(pst_2, dat = data)
trend_stacked(pst_3, dat = data_miss)

bind_rows(mutate(pst_1, method = 1),
          mutate(pst_2, method = 2),
          mutate(pst_3, method = 3)) %>%
  mutate(method = recode_factor(method,
                                `1` = "Original Method on Full Data",
                                `2` = "Multinomial Method on Full Data",
                                `3` = "Multinomial Method on Censored Data")) %>%
  group_by(method, iter) %>%
  summarize(llobs = unique(loglik_obs)/1500) %>%
  ggplot(aes(x = iter, y = llobs)) +
  facet_wrap(~method, ncol = 1, scales  = "free_y") +
  geom_line() +
  geom_point() +
  labs(y = "Observed Log Likelihood per observation") +
  theme_clean()



# plot max of diff
ggsave("figures/sim_EM_change-in-params_iia-method_missing.pdf", w = 4, h = 4)
ggsave("figures/sim_EM_change-in-params_iia-method_no-missing.pdf", w = 4, h = 4)



params <- read_rds("data/sim-data_missing.Rds")
params$theta

est <- store_iter[[40]]
round(est$theta, 2)

store_iia <- read_rds("data/EM/sim-iterations_IIA.Rds")
store_full <- read_rds("data/EM/sim-iterations.Rds")

for (t in length(store_iia)) {
  ord_iia <- order(store_iia[[t]]$theta, decreasing = TRUE)
  store_iia[[t]]$theta <- store_iia[[t]]$theta[ord_iia]
  store_iia[[t]]$zeta <-  store_iia[[t]]$zeta[, ord_iia]
  store_iia[[t]]$mu <-    store_iia[[t]]$mu[, ord_iia, ]

  ord_full <- order(store_full[[t]]$theta, decreasing = TRUE)
  store_full[[t]]$theta <- store_full[[t]]$theta[ord_full]
  store_full[[t]]$zeta <-  store_full[[t]]$zeta[, ord_full]
  store_full[[t]]$mu <-    store_full[[t]]$mu[ord_full, , ]
}

params$mu[1:3, 1, 3]
round(store_full[[40]]$mu[1:3, 1, 3], 2)
round(store_iia[[40]]$mu[1:3, 1, 3], 2)


cbind(params$Z, round(store_full[[40]]$zeta[, ], 2))[1:5, ]

tibble(as.numeric(params$theta), store_full[[40]]$theta)
round(store_iia[[40]]$zeta[1:5, ], 2)


