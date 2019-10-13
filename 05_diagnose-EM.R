library(tidyverse)
library(foreach)
library(ggthemes)
library(lemon)
library(glue)

source("03_define_EM-fun.R")

# Data --------
data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")
store_iter <- read_rds("data/EM/sim-iterations.Rds")


# Diagnostics ---------
params_stacked <- summ_params(store_iter, data = data, calc_loglik = FALSE)
params_stacked <- summ_params(store_iter, data = data, calc_loglik = TRUE) %>%
  mutate(llobs_scale = loglik_obs / (data$N*data$D))



# plot max of diff
params_stacked %>%
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
ggsave("figures/sim_EM_change-in-params_iia-method_no-missing.pdf", w = 4, h = 4)
