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
  L = n_distinct(unlist(unique_y_mat)) - 1,
  y = ch_sub[, ind_y],
  uy = unique_y_mat,
  n_u = n_u,
  U = U
)


# EM -----
source("03_define_EM-fun.R")
store_iter <- cat_mixture(data_ch, user_K = 3, n_iter = 50, fast = FALSE, IIA = TRUE)


# Store -------
write_rds(store_iter, "data/EM/charleston-iterations_trichotomous_iia.Rds")

#
# # Diagnose ----
# ch_iter <- read_rds("data/EM/charleston-iterations_trichotomous.Rds")
# ch_params <- summ_params(ch_iter, data_ch, calc_loglik = FALSE)
#
#
# ch_params %>%
#   mutate(llobs_scale = loglik_obs / (data_ch$N*data_ch$D)) %>%
#   group_by(iter) %>%
#   summarize(`Maximum Change in Parameter (probability scale)` = max(diff),
#             `Observed Log Likelihood (per data point)` = unique(llobs_scale)) %>%
#   pivot_longer(cols = -c(iter), names_to = "metric", values_to = "value") %>%
#   ggplot(aes(iter, value)) +
#   facet_rep_wrap(~metric, scales = "free_y", ncol = 1) +
#   coord_capped_cart(bottom='both', left = 'both') +
#   geom_point(size = 0.5) +
#   geom_line() +
#   theme_clean() +
#   theme(plot.background = element_rect(color = NA),
#         axis.line = element_line(color = "black"),
#         plot.title = element_text(hjust = 0.5, face = "bold"),
#         plot.caption = element_text(size = 6),
#         strip.background = element_rect(fill = "lightgray"))
# ggsave("figures/ch-complete_EM_change-in-params.pdf", w = 5, h = 3.5)
#
