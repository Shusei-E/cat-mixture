library(tidyverse)
library(klaR)
select <- dplyr::select

ch_split <- read_rds("data/split_ch-2018.Rds")


k_vanilla <- kmeans(na.omit(ch_split)[, 7:10], 5)

kdf_vanilla <- as_tibble(k_vanilla$centers) %>%
  mutate(n = k_vanilla$size)
kdf_vanilla


k_modes <- kmodes(na.omit(ch_split)[, 7:10], 5)
kdf_modes <- as_tibble(k_modes$modes) %>%
  mutate(n = as.integer(k_modes$size),
         cluster = 1:n()) %>%
  select(cluster, n, everything()) %>%
  arrange(desc(n))

kdf_modes

