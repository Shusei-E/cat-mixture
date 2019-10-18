library(tidyverse)
library(scales)
library(lemon)

# Data ------
ch_iter <- read_rds("data/EM/charleston-iterations_trichotomous.Rds")
ch_params <- summ_params(ch_iter[99:100], data_ch, calc_loglik = FALSE)


# vis params ---
ch_final <- summ_params(ch_iter, data_ch, calc_loglik = FALSE) %>% filter(iter == max(iter))

thetas <- ch_final %>%
  filter(str_detect(param_id, "theta")) %>%
  pull(values_now)

theta_df <- tibble(k = as.character(1:length(thetas)),
                   theta = thetas)

param_df <- ch_final %>%
  filter(str_detect(param_id, "mu")) %>%
  select(index = param_id, mu = values_now) %>%
  separate(index, into = c("param", "k", "j", "l"), sep = "_") %>%
  left_join(theta_df, by = "k") %>%
  mutate(k = fct_reorder(k, theta, .desc = TRUE),
         j = recode_factor(j,
                           `1` = "US House",
                           `2` = "St House",
                           `4` = "Probate Judge",
                           `3` = "County Council"),
         l = recode_factor(l,
                           `2` = "Straight",
                           `1` = "Split",
                           `0` = "Abstain")) %>%
  mutate(type = as.integer(k),
         type_lbl = glue("Type {type} ({percent(theta, accuracy = 1)})"))


param_df %>%
  ggplot() +
  geom_col(aes(x = fct_rev(j), y = mu, fill = fct_rev(l))) +
  scale_fill_viridis_d(name = NULL) +
  coord_capped_flip(bottom = "none",
                    left = brackets_vertical(direction = "right")) +
  facet_rep_wrap(~type_lbl) +
  theme_light() +
  labs(x = "", y = "Probability of  voting ...") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.position = "bottom")

ggsave("figures/ch_types-described.pdf", w = 7, h = 2.5)
