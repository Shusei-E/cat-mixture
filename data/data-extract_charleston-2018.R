library(tidyverse)
library(haven)

wide <- read_dta("data/output/by-person_votechoice.dta")

ch_wide <- wide %>% 
  filter(county == "Charleston",
         elec == "2018-11-06")

ch_sel <- ch_wide %>% 
  filter(USH_dist == 1) %>% 
  select_if(function(x) sum(!is.na(x)) > 0)


write_rds(ch_sel, "sims/data/charleston-2018.Rds")
