library(tidyverse)

ch_raw <- read_rds("data/charleston-2018.Rds")


ch_votes <- ch_raw %>%
  select(1:5, matches("party")) %>%
  select(-PTY_party)


ch_split <- ch_votes %>%
  rename_at(vars(matches("(USH|HOU|CCD|JPR)_party")), ~str_replace(.x, "_party", "_split")) %>%
  mutate_at(
    vars(matches("_split")),
    ~case_when(GOV_party == .x & GOV_party %in% c(-1, 1) & .x %in% c(-1, 1) ~ 2, # straight
               GOV_party != .x & GOV_party %in% c(-1, 1)  & .x %in% c(-1, 1) ~ 1, # split
               .x == 0 & GOV_party %in% c(-1, 1) ~ 0, # abstain
               is.na(GOV_party) | GOV_party == 0 ~ NA_real_,
               is.na(.x) ~ NA_real_,
               TRUE ~ 0)
    ) %>%
  mutate(GOV_party = recode(GOV_party, `1` = "R", `-1` = "D", `0` = "A"))


write_rds(ch_split, "data/split_ch-2018.Rds")
