library(tidyverse)
library(patchwork)
library(ggthemes)
library(lemon)
library(glue)

source("03_define_EM-fun.R")

data <- read_rds("data/sim-data.Rds")
data_miss <- read_rds("data/sim-data_missing.Rds")

# iters
fit_k1 <- read_rds("data/EM/n300_kmeans-init/sim-iter.Rds")
fit_k3 <- read_rds("data/EM/n300_kmeans-init/sim-iter_IIA-miss.Rds")

fit_e1 <- read_rds("data/EM/sim-iterations.Rds")
fit_e3 <- read_rds("data/EM/sim-iterations_IIA-miss.Rds")

fit_e3[[1]]$opts

# estimated params
pst_1 <- read_rds("data/EM/sim-stats.Rds")
pst_3 <- read_rds("data/EM/sim-stats_IIA-miss.Rds")

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