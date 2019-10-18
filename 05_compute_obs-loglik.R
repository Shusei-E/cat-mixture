library(tidyverse)
library(foreach)

source("03_define_EM-fun.R")

# Data --------
data <- read_rds("data/sim-data.Rds")
data_miss <- read_rds("data/sim-data_missing.Rds")
params <- read_rds("data/sim-params.Rds")
store_iter_1 <- read_rds("data/EM/sim-iterations.Rds")
store_iter_2 <- read_rds("data/EM/sim-iterations_IIA-full.Rds")
store_iter_3 <- read_rds("data/EM/sim-iterations_IIA-miss.Rds")


# Calculate metrics  ---------
pst_3_new <- summ_params(store_iter_3, data = data_miss, IIA = TRUE, calc_loglik = TRUE)
write_rds(pst_3_new, "data/EM/sim-stats_IIA-miss.Rds")

pst_1_new <- summ_params(store_iter_1, data = data, calc_loglik = TRUE)
write_rds(pst_1_new, "data/EM/sim-stats.Rds")
pst_2_new <- summ_params(store_iter_2, data = data, IIA = TRUE, calc_loglik = TRUE)
write_rds(pst_2_new, "data/EM/sim-stats_IIA-full.Rds")