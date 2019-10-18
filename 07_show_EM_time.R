store_iter_3 <- read_rds("data/EM/sim-iterations_IIA-miss.Rds")
pst <- read_rds("data/EM/sim-stats_IIA-miss.Rds")

pst %>%
  group_by(iter) %>%
  ggplot(aes(x = iter, y = time_l)) +
  geom_point() +
  geom_line()

pst_L <- pst %>%
  distinct(iter, L = time_l)

store_iter_3 %>%
  map_dfr(~tibble(E = .x$opts$time_e, M = .x$opts$time_m)) %>%
  mutate(iter = 1L:n()) %>%
  left_join(pst_L, by = "iter") %>%
  # summarize_all(~mean(.x, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(iter, E)) +  geom_line()



# n = 3000, D = 5: 300 seconds per iteration
# n = 300, D= 5: 28-30 seconds per iteration
#