rm(list = ls())

library(tidyverse)
library(RWiener)
library(rutils)
library(grid)
library(gridExtra)

tbl_transfer <- readRDS("data/transfer-rt-data.RDS")
tbl_transfer$label <- "new"
tbl_old <- readRDS("data/hotspot-data.RDS")
tbl_old$label <- "old"

tbl_test <- rbind(
  tbl_transfer %>% dplyr::select(x1, x2, label), 
  tbl_old %>% dplyr::select(x1, x2, label)
)

params_fin <- c(1.23, .5, .5)

sims_all <- pmap_dbl(
  tbl_test[, c("x1", "x2")], 
  ~ sum(pmap_dbl(
    tbl_old[, c("x1", "x2")], 
    f_similarity, 
    c(params_fin[[2]], 1 - params_fin[[2]]), 15, tibble(x1 = .x, x2 = .y), 1
  ))
)
p_gamma <- .25
p_thx_response <- .2
sims_all_resp_prop <- sims_all^p_gamma / max(sims_all^p_gamma)


tbl_test$sim_all <- sims_all
tbl_test$sim_all_resp_prop <- sims_all_resp_prop


ggplot(tbl_test, aes(sim_all_resp_prop)) +
  geom_histogram(fill = "white", color = "red") +
  facet_wrap(~ label, scales = "free_y")



alpha <- 1.23
beta <- .55
delta_ic <- 2#.93
delta_slope <- -.02
tau <- .27

v_delta_slopes <- seq(-1, 1, by = .2)

l_out <- list()
n_reps <- 100

for (i in 1:length(v_delta_slopes)) {
  
  tbl_rt_strat <- map_df(
    scale(tbl_transfer$sim_strat, center = TRUE, scale = FALSE), 
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
      sim_strat = rep(tbl_transfer$sim_strat, each = n_reps),
      sim_hotspot = rep(tbl_transfer$sim_hotspots, each = n_reps)
    )
  tbl_rt_strat$item_id <- rep(1:nrow(tbl_transfer), each = n_reps)
  
  tbl_rt_hotspots <- map_df(
    scale(tbl_transfer$sim_hotspots, center = TRUE, scale = FALSE), 
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
      sim_strat = rep(tbl_transfer$sim_strat, each = n_reps),
      sim_hotspot = rep(tbl_transfer$sim_hotspots, each = n_reps)
    )
  tbl_rt_hotspots$item_id <- rep(1:nrow(tbl_transfer), each = n_reps)
  
  tbl_both <- tbl_rt_strat %>% dplyr::select(sim_strat, rt_strat = rt, resp_strat = resp, item_id) %>%
    cbind(
      tbl_rt_hotspots %>% 
        dplyr::select(sim_hotspot, rt_hotspot = rt, resp_hotspot = resp)
    ) %>%
    mutate(
      resp_same = resp_strat == resp_hotspot,
      rt_diff = rt_strat - rt_hotspot
    )
  
  tbl_both_agg <- grouped_agg(tbl_both, c(item_id, sim_strat, sim_hotspot), c(resp_same, rt_diff))
  
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



ggplot(tbl_both_agg, aes(sim_strat - sim_hotspot, mean_rt_diff)) +
  geom_line()



ggplot(tbl_both_agg, aes(sim_strat - sim_hotspot, mean_resp_same)) +
  geom_line()

ggplot(tbl_both_agg, aes(sim_strat, sim_hotspot)) + geom_point() + geom_abline()
