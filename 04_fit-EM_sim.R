library(tidyverse)
library(foreach)
library(glue)

# Functions ---
source("03_define_EM-fun.R")

# No missing
data <- read_rds("data/sim-data.Rds")
store_iter <- cat_mixture(data, user_K = 3, n_iter = 50, fast = FALSE, IIA = TRUE)
write_rds(store_iter, "data/EM/sim-iterations.Rds")


data <- read_rds("data/sim-data_missing.Rds")
store_iter <- cat_mixture(data, user_K = 3, n_iter = 25, fast = FALSE, IIA = TRUE)
write_rds(store_iter, "data/EM/sim-iterations_IIA.Rds")

