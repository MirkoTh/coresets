rm(list = ls())
set.seed(3453)

library(tidyverse)
library(grid)
library(gridExtra)
library(furrr)
library(docstring)
library(rutils)
library(cmdstanr)
library(RWiener)

path_load <- c("utils/utils.R", "scripts/stan-wiener.R")
walk(path_load, source)


fl_pths <- str_c("data/cd-cr/", dir("data/cd-cr/"))
fl_pths <- fl_pths[fl_pths != "data/cd-cr/Data_S4_4.txt"]
l_tbl_raw <- map(fl_pths[endsWith(fl_pths, ".txt")], read.table, header = TRUE)
tbl_raw <- reduce(l_tbl_raw, rbind) %>% as_tibble()
names(tbl_raw)[1:15] <- c("subj", "seed", "sess", "block", "trial", "taskType", "testType", "resp", "rt",
                          "setSize", "targLoc", "targValue", "recogValue", "respValue", "deviance")

tbl_raw$taskType <- as.factor(tbl_raw$taskType)
levels(tbl_raw$taskType) <- c("Recognition","Recall")
tbl_raw$c_size <- tbl_raw$targValue-tbl_raw$recogValue
tbl_raw$c_size[tbl_raw$c_size >= 180] <- -(360 - tbl_raw$c_size[tbl_raw$c_size >= 180])
tbl_raw$c_size[tbl_raw$c_size < -180] <- -(-360 - tbl_raw$c_size[tbl_raw$c_size < -180])
tbl_raw$c_size[tbl_raw$taskType == "Recall"] <- -999
tbl_raw$c_size_abs <- abs(tbl_raw$c_size)


# trimming
use2 <- tbl_raw$rt > 0.18 & tbl_raw$rt < 2.5
use1 <- rep(0, length(tbl_raw[, 11]))
ctrial <- 0

for (subj in 1:max(tbl_raw$subj)){
  tmp <- tbl_raw$rt[tbl_raw[,1] == subj & use2]; mn <- mean(tmp); sds <- sd(tmp); temp <- mn+2.5*sds	
  use1[tbl_raw[, 1] == subj & tbl_raw$rt < temp] <- 1
  ctrial <- c(ctrial,1:nrow(tbl_raw[tbl_raw[, 1] == subj, 1]))
}
ctrial <- ctrial[-1]
use <- use1 & use2 & ctrial > 50
tbl_raw <- tbl_raw[use, ]

# for this initial wiener analysis, exclude subj 3
ggplot(tbl_raw %>% filter(taskType == "Recognition"), aes(rt)) +
  geom_histogram() +
  facet_wrap(~ subj)

tbl_cd <- tbl_raw %>% filter(taskType == "Recognition" & subj != 3 & rt > .3)
tbl_cr <- tbl_raw %>% filter(taskType == "Recall" & subj != 3)

tbl_cd$id_rand <- sample(1:nrow(tbl_cd), nrow(tbl_cd), replace = FALSE)
tbl_cd$resp_recode <- map_chr(tbl_cd$resp + 1, ~ c("upper", "lower")[.x])

tbl_cd %>% group_by(c_size_abs = abs(c_size)) %>% summarize(acc = mean(resp), m_rt = mean(rt), n = n()) %>%
  ggplot(aes(c_size_abs, acc)) +
  geom_line() +
  geom_point(aes(size = n))


# 
# l_data <- list(
#   n_data_0 = nrow(tbl_cd %>% filter(resp == 0)),
#   rt_0 = tbl_cd$rt[tbl_cd$resp == 0] * 1,
#   n_data_1 = nrow(tbl_cd %>% filter(resp == 1)),
#   rt_1 = tbl_cd$rt[tbl_cd$resp == 1] * 1
# )
# 
# txt_wiener_base <- stan_wiener()
# m_wiener_base <- cmdstan_model(txt_wiener_base)
# 
# init_fun <- function() list(delta = .01)
# 
# 
# fit_wiener_base <- m_wiener_base$sample(
#   data = l_data, iter_sampling = 200, iter_warmup = 100, chains = 1#, init = init_fun
# )
# 
# 
# 
# pars_interest <- c("alpha", "beta", "delta", "tau")
# tbl_draws <- fit_wiener_base$draws(variables = pars_interest, format = "df")
# tbl_summary <- fit_wiener_base$summary()#variables = pars_interest)
# 
# tbl_posterior <- tbl_draws %>% 
#   as_tibble() %>%
#   dplyr::select(c(all_of(pars_interest), ".chain")) %>%
#   rename(chain = .chain) %>%
#   pivot_longer(-chain, names_to = "parameter", values_to = "value") %>%
#   mutate(parameter = factor(parameter))
# 
# ggplot(tbl_posterior, aes(value)) +
#   geom_histogram() +
#   facet_wrap(~ parameter, scales = "free_x")
# 
# 
# params_bf <- c("Intercept", "Trial (Binned)")
# l <- sd_bfs(tbl_posterior, params_bf, sqrt(2)/4)
# bfs <- l[[1]]
# tbl_thx <- l[[2]]
# 
# # plot the posteriors and the bfs
# map(as.list(params_bf), plot_posterior, tbl_posterior, tbl_thx, bfs)
# 

# RWiener
# 
# set.seed(0)
# dat <- rwiener(n=100, alpha=2, tau=.3, beta=.5, delta=.5)
# optim1 <- optim(c(1, .1, .1, 1), wiener_deviance, dat=dat, method="Nelder-Mead")
# 
# dat <- rwiener(n=10000, alpha=2, tau=.3, beta=.5, delta=.5)
# optim2 <- optim(c(1, .1, .1, 1), wiener_deviance, dat=dat, method="Nelder-Mead")
# 
# 
# dwiener(.75, 2, .3, .5, .5)
# plot(dwiener(seq(0, 5, by = .01), 2, .3, .5, .5, resp = rep("upper", 501)))
# 


wiener_reg_delta_log <- function(x, my_tbl) {
  #' @description Wiener LL with linear regression on drift rate
  
  alpha <- x[[1]]
  tau <- x[[2]]
  beta <- x[[3]]
  delta_ic <- x[[4]]
  delta_slope <- x[[5]]
  
  lik <- pmap_dbl(
    my_tbl[, c("rt", "resp_recode", "pred_lr")], ~ dwiener(
      q = ..1, alpha = alpha, tau = tau, beta = beta, 
      delta = delta_ic + delta_slope * ..3, resp = ..2
    )
  )
  
  neg2loglik <- -2*sum(log(pmax(lik,1e-10)))
  
  return(neg2loglik)
}


tbl_cd %>% filter(resp == 1 & c_size > 0) %>%
  mutate(c_size_cut = cut(c_size, 7)) %>%
  group_by(c_size_cut) %>%
  summarize(m_rt = mean(rt), n = n()) %>%
  ggplot(aes(c_size_cut, m_rt, group = 1)) + 
  geom_line() + 
  geom_point(aes(size = n))


results_c_size <- optim(c(1, .1, .1, 1, .05), wiener_reg_delta_log, my_tbl = tbl_cd[1:1000, ] %>% mutate(pred_lr = c_size_abs))


# now use these parameters to predict rt data and choices for the same stimuli used in the two-dimensional categorization task
rwiener(1, results_c_size$par[[1]], results_c_size$par[[2]], results_c_size$par[[3]], results_c_size$par[[4]] + results_c_size$par[[5]] * c(1, 2))
l_tbl_up_and_down <- readRDS("data/l_tbl_up_and_down.RDS")
tbl_hotspots <- readRDS("data/hotspot-data.RDS")
tbl_transfer <- readRDS("data/transfer-data.RDS")

params_fin <- c(1.23, .5, .5)
lo <- c(0, .01, .0001, 1e-10)
hi <- c(10, .99, .9999, .999999999)

params <- pmap(list(params_fin[1:3], lo[1:3], hi[1:3]), upper_and_lower_bounds)

sims_hotspots <- pmap_dbl(
  tbl_transfer[, c("x1", "x2")], 
  ~ sum(pmap_dbl(
    tbl_hotspots[, c("x1", "x2")], 
    f_similarity, 
    c(params_fin[[2]], 1 - params_fin[[2]]), params_fin[[3]], tibble(x1 = .x, x2 = .y), 1
  ))
)

sims_strat <- pmap_dbl(
  tbl_transfer[, c("x1", "x2")], 
  ~ sum(pmap_dbl(
    l_tbl_up_and_down[[7]][, c("x1", "x2")], 
    f_similarity, 
    c(params_fin[[2]], 1 - params_fin[[2]]), params_fin[[3]], tibble(x1 = .x, x2 = .y), 1
  ))
)

p_gamma <- 3
p_thx_response <- .2
sims_strat <- sims_strat^p_gamma / max(sims_strat^p_gamma)
sims_hotspots <- sims_hotspots^p_gamma / max(sims_hotspots^p_gamma)


plot(sort(sims_strat)/(sort(sims_strat) + p_thx_response))
plot(sort(sims_hotspots)/(sort(sims_hotspots) + p_thx_response))


tbl_transfer$sim_strat <- sims_strat
tbl_transfer$sim_hotspots <- sims_hotspots
tbl_transfer$pred_diff <- tbl_transfer$sim_strat - tbl_transfer$sim_hotspots
tbl_transfer$pred_diff_abs <- abs(tbl_transfer$pred_diff)
tbl_transfer$disc_power <- rank(desc(tbl_transfer$pred_diff_abs))

tbl_transfer_disc <- tbl_transfer %>% filter(disc_power <= 10)

ggplot(tbl_transfer_disc, aes(x1, x2)) +
  geom_point(aes(size = pred_diff, color = pred_diff)) + 
  ggrepel::geom_label_repel(aes(label = round(pred_diff, 2))) +
  theme_bw() +
  scale_x_continuous(expand = c(.03, 0)) +
  scale_y_continuous(expand = c(.03, 0)) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme(strip.background = element_rect(fill = "white")) +
  scale_size_continuous(name = "Pred. Difference", range = c(1, 5)) +
  scale_color_gradient2(midpoint = 0, low = "skyblue2", high = "tomato4", guide = "none") +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 7))

n_reps <- 20

tbl_rt_strat <- map_df(
  tbl_transfer$sim_strat, 
  ~ rwiener(
    n = n_reps, alpha = results_c_size$par[[1]], tau = results_c_size$par[[2]], 
    beta = results_c_size$par[[3]], delta = results_c_size$par[[4]] + .05*.x#results_c_size$par[[5]] * .x
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

tbl_rt_hotspots <- map_df(
  tbl_transfer$sim_hotspots, 
  ~ rwiener(
    n = n_reps, alpha = results_c_size$par[[1]], tau = results_c_size$par[[2]], 
    beta = results_c_size$par[[3]], delta = results_c_size$par[[4]] + .05*.x#results_c_size$par[[5]] * .x
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

tbl_rt_both <- rbind(tbl_rt_strat, tbl_rt_hotspots)

tbl_acc <- grouped_agg(tbl_rt_both %>% mutate(p_old = resp == "old"), c(model, resp), c(rt)) %>%
  mutate(prop_resp = n / sum(n))

ggplot(rbind(tbl_rt_hotspots, tbl_rt_strat), aes(rt)) +
  geom_histogram(fill = "skyblue2", color = "black") +
  geom_label(data = tbl_acc, aes(x = 1.5, y = 40, label = str_c("Prop. = ", round(prop_resp, 3)))) +
  facet_grid(model ~ resp) + 
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "RT (s)", y = "Nr. Responses") +
  theme(strip.background = element_rect(fill = "white"))

` 
++`tbl_ll_diff <- tibble(n_reps = numeric(), it = numeric(), ll_diff = numeric())


for (nr in c(5, 15, 30)) {
  for (i in 1:50) {
    tbl_rt_strat <- map_df(
      tbl_transfer_disc$sim_strat, 
      ~ rwiener(
        n = nr, alpha = results_c_size$par[[1]], tau = results_c_size$par[[2]], 
        beta = results_c_size$par[[3]], delta = results_c_size$par[[4]] + 3 * .x
      )
    ) %>% as_tibble() %>% 
      rename(rt = q) %>%
      mutate(
        model = "Strat. Sampling",
        resp_recode = resp,
        resp = fct_relabel(resp, ~ c("old", "new")),
        sim_strat = rep(tbl_transfer_disc$sim_strat, each = nr),
        sim_hotspot = rep(tbl_transfer_disc$sim_hotspots, each = nr)
      )
    
    results_recover_strat_by_strat <- optim(
      c(1, .1, .1, 1, .05), 
      wiener_reg_delta_log, 
      my_tbl = tbl_rt_strat %>% mutate(pred_lr = sim_strat)
    )
    results_recover_strat_by_hotspot <- optim(
      c(1, .1, .1, 1, .05), 
      wiener_reg_delta_log, 
      my_tbl = tbl_rt_strat %>% mutate(pred_lr = sim_hotspot)
    )
    tbl_ll_diff <- rbind(
      tbl_ll_diff,
      tibble(
        n_reps = nr,
        it = i,
        ll_diff = results_recover_strat_by_strat$value - results_recover_strat_by_hotspot$value
      )
    )
    saveRDS(tbl_ll_diff, file = "data/wiener-strat-sampling-recovery-disc-power-top10.RDS")
    
  }
}

tbl_ll_diff <- readRDS("data/wiener-strat-sampling-recovery-disc-power-top10.RDS")


ggplot(tbl_ll_diff, aes(ll_diff)) +
  geom_histogram(color = "black", fill = "skyblue2") +
  geom_vline(xintercept = 0) +
  facet_wrap(~ n_reps) + 
  theme_bw() +
  scale_x_continuous(expand = c(.1, 0)) +
  scale_y_continuous(expand = c(.1, 0)) +
  labs(x = "LL Strat. Sampling - LL Hotspots", y = "Nr. Samples", title = "Strat. Sampling Generating") +
  theme(strip.background = element_rect(fill = "white"))
  


# todos
# use a mixture from training set (actually old responses) and transfer set (new responses)
# generate this mixture using data points from the transfer set, on which predictions from the two models differ a lot




