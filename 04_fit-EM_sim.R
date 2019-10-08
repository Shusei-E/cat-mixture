library(tidyverse)
library(foreach)
library(glue)

# Functions ---
source("03_define_EM-fun.R")

# Data --------
data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")

# store sim
store_iter <- cat_mixture(data, user_K = 3, n_iter = 3)
write_rds(store_iter, "data/EM/sim-iterations.Rds")