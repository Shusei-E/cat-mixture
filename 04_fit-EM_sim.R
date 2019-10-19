library(tidyverse)
library(foreach)
library(glue)

# Functions ---
source("03_define_EM-fun.R")


# default
data <- read_rds("data/sim-data.Rds")
store_iter_rrr <- cat_mixture(data, user_K = 3, n_iter = 3, fast = TRUE, init = "equal")
store_iter_cpp <- cat_mixture(data, user_K = 3, n_iter = 3, fast = TRUE, init = "equal")
write_rds(store_iter, "data/EM/n300_kmeans-init/sim-iterations.Rds")

store_iter_rrr[[2]]$opts[[6]]
store_iter_cpp[[2]]$opts[[6]]

as.numeric(store_iter_rrr[[2]]$opts[[6]]) /
as.numeric(store_iter_cpp[[2]]$opts[[6]])

store_iter_rrr[[3]]$mu[, , 2]
store_iter_cpp[[3]]$mu[, , 2]

# using mlogit algorithm on missing data
data_miss <- read_rds("data/sim-data_missing.Rds")
store_iter <- cat_mixture(data_miss, user_K = 3, n_iter = 50, IIA = TRUE, init = "kmeans")
write_rds(store_iter, "data/EM/n300_kmeans-init/sim-iterations_IIA-miss.Rds")


# using mnlogit algorithm on full data
# store_iter <- cat_mixture(data, user_K = 3, n_iter = 50, IIA = TRUE)
# write_rds(store_iter, "data/EM/sim-iterations_IIA-full.Rds")

