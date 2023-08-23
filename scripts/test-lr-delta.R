rm(list = ls())


library(tidyverse)
library(RWiener)
library(rutils)
library(grid)
library(gridExtra)


dirs_home_grown <- c("utils/utils.R")
walk(dirs_home_grown, source)


tbl_transfer <- readRDS("data/transfer-rt-data.RDS")
tbl_transfer$label <- "new"
tbl_old <- readRDS("data/hotspot-data.RDS")
tbl_old$label <- "old"

l_tbl_strat <- readRDS("data/l_tbl_up_and_down.RDS")
tbl_strat <- l_tbl_strat[[2]]

tbl_test <- rbind(
  tbl_transfer %>% dplyr::select(x1, x2, label), 
  tbl_old %>% dplyr::select(x1, x2, label)
)



params_fin <- c(1.23, .5, .5)

sims_hotspot <- pmap_dbl(
  tbl_test[, c("x1", "x2")], 
  ~ sum(pmap_dbl(
    tbl_old[, c("x1", "x2")], 
    f_similarity, 
    c(params_fin[[2]], 1 - params_fin[[2]]), 15, tibble(x1 = .x, x2 = .y), 1
  ))
)

sims_strat <- pmap_dbl(
  tbl_test[, c("x1", "x2")], 
  ~ sum(pmap_dbl(
    tbl_strat[, c("x1", "x2")], 
    f_similarity, 
    c(params_fin[[2]], 1 - params_fin[[2]]), 15, tibble(x1 = .x, x2 = .y), 1
  ))
)

cor(scale(sims_hotspot)[, 1], scale(sims_strat)[, 1])
plot(scale(sims_hotspot)[, 1], scale(sims_strat)[, 1])
hist(sims_strat)
hist(sims_hotspot)

p_gamma <- 1
p_thx_response <- .2

# the two similarities are on different scales because of different number of exemplars in memory
# do not bother about the absolute values of the two similarities 
# because they are inserted as regressors onto drift rate
# by scaling both similarities, we can compare the coefficients how strong they affect v

sims_hotspot_resp_prop <- scale(sims_hotspot^p_gamma)[, 1] #/ max(sims_hotspot^p_gamma)
sims_strat_resp_prop <- scale(sims_strat^p_gamma)[, 1] #/ max(sims_strat^p_gamma)

cor(sims_hotspot_resp_prop, sims_strat_resp_prop)
plot(sims_hotspot_resp_prop, sims_strat_resp_prop)

tbl_test$sim_hotspots <- sims_hotspot
tbl_test$sim_hotspots_resp_prop <- sims_hotspot_resp_prop
tbl_test$sim_strat <- sims_strat
tbl_test$sim_strat_resp_prop <- sims_strat_resp_prop


ggplot(tbl_test %>% pivot_longer(c(sim_hotspots, sim_strat)), aes(value)) +
  geom_histogram(fill = "white", color = "red") +
  facet_grid(name ~ label, scales = "free")


ggplot(tbl_test %>% pivot_longer(c(sim_hotspots_resp_prop, sim_strat_resp_prop)), aes(value)) +
  geom_histogram(fill = "white", color = "red") +
  facet_grid(name ~ label, scales = "free")


alpha <- 1.23
beta <- .55
delta_ic <- 2#.93
delta_slope <- -.02
tau <- .27

v_delta_slopes <- seq(-3, 3, by = .6)

l_out <- list()
n_reps <- 100

for (i in 1:length(v_delta_slopes)) {
  
  tbl_rt_strat <- map_df(
    scale(tbl_test$sim_strat_resp_prop, center = TRUE, scale = FALSE), 
    ~ rwiener(
      n = n_reps, alpha = alpha, tau = tau, beta = beta,
      delta = delta_ic + v_delta_slopes[[i]] * .x #results_c_size$par[[5]] * .x
    )
  ) %>% as_tibble() %>% 
    rename(rt = q) %>%
    mutate(
      model = "Strat. Sampling",
      resp_recode = resp,
      resp = fct_relabel(resp, ~ c("old", "new")),
      sim_strat = rep(tbl_test$sim_strat_resp_prop, each = n_reps),
      sim_hotspots = rep(tbl_test$sim_hotspots_resp_prop, each = n_reps)
    )
  tbl_rt_strat$item_id <- rep(1:nrow(tbl_test), each = n_reps)
  
  tbl_rt_hotspots <- map_df(
    scale(tbl_test$sim_hotspots_resp_prop, center = TRUE, scale = FALSE), 
    ~ rwiener(
      n = n_reps, alpha = alpha, tau = tau, beta = beta,
      delta = delta_ic + v_delta_slopes[[i]] * .x #results_c_size$par[[5]] * .x
    )
  ) %>% as_tibble() %>% 
    rename(rt = q) %>%
    mutate(
      model = "Hotspot",
      resp_recode = resp,
      resp = fct_relabel(resp, ~ c("old", "new")),
      sim_strat = rep(tbl_test$sim_strat_resp_prop, each = n_reps),
      sim_hotspots = rep(tbl_test$sim_hotspots_resp_prop, each = n_reps)
    )
  tbl_rt_hotspots$item_id <- rep(1:nrow(tbl_test), each = n_reps)
  
  tbl_both <- tbl_rt_strat %>% 
    dplyr::select(sim_strat, rt_strat = rt, resp_strat = resp, item_id) %>%
    cbind(
      tbl_rt_hotspots %>% 
        dplyr::select(sim_hotspots, rt_hotspot = rt, resp_hotspot = resp)
    ) %>%
    mutate(
      resp_same = resp_strat == resp_hotspot,
      rt_diff = rt_strat - rt_hotspot
    )
  
  tbl_both_agg <- grouped_agg(
    tbl_both, 
    c(item_id, sim_strat, sim_hotspots), 
    c(resp_same, rt_diff, rt_strat)
  )
  
  tbl_rt_both <- rbind(tbl_rt_strat, tbl_rt_hotspots)
  
  tbl_acc <- grouped_agg(tbl_rt_both %>% mutate(p_old = resp == "old"), c(model, resp), c(rt)) %>%
    mutate(prop_resp = n / sum(n)) %>% ungroup()
  
  tbl_eval <- pivot_wider(tbl_acc, id_cols = resp, names_from = model, values_from = prop_resp) %>%
    mutate(var = "prop_resp") %>%
    rbind(
      pivot_wider(tbl_acc, id_cols = resp, names_from = model, values_from = mean_rt) %>%
        mutate(var = "rt")
    )
  
  tbl_eval <- tbl_eval %>% mutate(diff = Hotspot - `Strat. Sampling`)
  
  diff_prop_resp <- tbl_eval %>% filter(resp == "old" & var == "prop_resp") %>%
    dplyr::select(diff) %>% as_vector()
  
  diff_rts <- tbl_eval %>% filter(var == "rt") %>%
    dplyr::select(diff) %>% as_vector()
  
  l_out[[i]] <- tbl_both_agg
  
}


tbl_abs <- as.data.frame(reduce(map(l_diffs, "Hotspot"), rbind)) %>% 
  mutate(model_name = "Hotspot") %>%
  rbind(
    as.data.frame(reduce(map(l_diffs, "StrSam"), rbind)) %>% 
      mutate(model_name = "Strat. Sampl.") 
  ) %>%
  as_tibble() %>%
  mutate(delta_slope = rep(v_delta_slopes, 2))
colnames(tbl_abs)[1:2] <- c("prop_old", "rt_old")


tbl_results <- tibble(
  slp = v_delta_slopes,
  diff_resp = map_dbl(l_diffs, "resp"),
  diff_rt_old = map_dbl(map(l_diffs, "rt"), 1),
  diff_rt_new = map_dbl(map(l_diffs, "rt"), 2)
)

tbl_results %>% pivot_longer(cols = c(diff_resp, diff_rt_old, diff_rt_new)) %>%
  ggplot(aes(slp, value, group = name)) +
  geom_line() +
  geom_point(size = 4, color = "white") +
  geom_point(aes(color = name)) +
  facet_wrap(~ name) +
  theme_bw() +
  scale_x_continuous(expand = c(.02, 0)) +
  scale_y_continuous(expand = c(.02, 0)) +
  labs(x = "Slope", y = "Prediction Difference") + 
  theme(strip.background = element_rect(fill = "white")) + 
  scale_color_viridis_d(name = "Variable")

ggplot(tbl_abs %>% pivot_longer(cols = c(prop_old, rt_old)), aes(delta_slope, value, group = name)) +
  geom_line(aes(color = name)) +
  geom_point(color = "white", size = 4) +
  geom_point(aes(color = name)) +
  theme_bw() +
  facet_wrap(~ model_name) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "") + 
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_manual(values = c("skyblue2", "tomato4"), name = "")



ggplot(tbl_both_agg, aes(sim_strat - sim_hotspots, mean_rt_diff)) +
  geom_line()
ggplot(tbl_both_agg, aes(sim_strat, mean_rt_strat)) +
  geom_line()



ggplot(tbl_both_agg, aes(sim_strat - sim_hotspots, mean_resp_same)) +
  geom_line()

ggplot(tbl_both_agg, aes(sim_strat, sim_hotspots)) + geom_point() + geom_abline()
