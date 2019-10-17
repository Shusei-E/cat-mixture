library(tidyverse)
library(foreach)
library(glue)

# Functions ---
source("03_define_EM-fun.R")



# using mlogit algorithm on missing data
data_miss <- read_rds("data/sim-data_missing.Rds")
store_iter <- cat_mixture(data_miss, user_K = 2, n_iter = 50, IIA = TRUE)
write_rds(store_iter, "data/EM/sim-iterations_IIA-miss.Rds")


# using mnlogit algorithm on full data
data <- read_rds("data/sim-data.Rds")
store_iter <- cat_mixture(data, user_K = 2, n_iter = 50, IIA = TRUE)
write_rds(store_iter, "data/EM/sim-iterations_IIA-full.Rds")

# default
store_iter <- cat_mixture(data, user_K = 4, n_iter = 50)
write_rds(store_iter, "data/EM/sim-iterations.Rds")
