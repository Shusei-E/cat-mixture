library(tidyverse)
library(haven)
library(fs)
lr <- "~/Dropbox/EL155s"
wide <- read_rds(path(lr, "output/04_wide/wide_coded.Rds"))

ch_wide <- wide %>%
  filter(county == "Charleston",
         elec == "2018-11-06")

ch_sel <- ch_wide %>%
  filter(USHOU_dist == 1) %>%
  select_if(function(x) sum(!is.na(x)) > 0) %>%
  select(-ctydist, -elec, -county, -ballot_style, -precinct_id) %>%
  rename_all(~str_replace_all(.x, "USHOU", "USH")) %>%
  rename_all(~str_replace_all(.x, "JPRB", "JPR"))

write_rds(ch_sel, "data/charleston-2018.Rds")
