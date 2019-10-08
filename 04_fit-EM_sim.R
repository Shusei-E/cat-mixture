library(tidyverse)
library(foreach)
library(glue)

# Functions ---
source("03_define_EM-fun.R")

# Data --------
data <- read_rds("data/sim-data.Rds")
params <- read_rds("data/sim-params.Rds")

# Unique profile and speed up? ----
fast <- TRUE
if (!fast) {
  data$U <- data$N
  data$n_u <- rep(1, data$N)
  data$uy <- data$y
}


store_iter <- cat_mixture(data, user_K = 3, n_iter = 3)
