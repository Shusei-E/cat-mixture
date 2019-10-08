library(tidyverse)

ch <- read_rds("data/split_ch-2018.Rds")

ch_sub <- ch %>% filter(!is.na(HOU_split),
              !is.na(CCD_split))

ind_y <- str_which(colnames(ch_sub), "_split")
D <- length(ind_y)

# unique profiles
profile_y <- ch_sub[, ind_y] %>%
  mutate(voter = 1:n()) %>%
  pivot_longer(-voter, names_to = "j", values_to = "y") %>%
  group_by(voter) %>%
  summarize(profile = str_c(y, collapse = "")) %>%
  count(profile, name = "n_u")

unique_y <- profile_y %>%
  separate(profile, into = str_c("v", 1:D), sep = 1:D) %>%
  mutate_all(as.integer)

unique_y_mat <- unique_y[, 1:D]
n_u <- unique_y$n_u
U <- length(n_u)


data_ch <- list(
  D = D,
  N = nrow(ch_sub),
  L = length(unique(unlist(unique_y_mat))),
  y = ch_sub[, ind_y],
  uy = unique_y_mat,
  n_u = n_u,
  U = U
)



# EM -----
source("03_define_EM-fun.R")
store_iter <- cat_mixture(data_ch, user_K = 2, n_iter = 100)

round(store_iter[[5]]$theta, 2)

# Store -------
write_rds(store_iter, "data/EM/charleston-iterations_trichotomous.Rds")

# ggsave("figures/ch-complete_EM_change-in-params.pdf", w = 4, h = 4)
