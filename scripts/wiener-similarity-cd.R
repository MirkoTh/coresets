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

tbl_cd %>% group_by(c_size_abs = abs(c_size)) %>% 
  summarize(p_new = mean(resp), m_rt = mean(rt), n = n()) %>%
  ggplot(aes(c_size_abs, p_new)) +
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


wiener_reg1_delta_log <- function(x, my_tbl) {
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

wiener_reg2_delta_log <- function(x, my_tbl) {
  #' @description Wiener LL with linear regression on drift rate
  
  alpha <- x[["alpha"]]
  tau <- x[["tau"]]
  beta <- x[["beta"]]
  delta_ic <- x[["delta_ic"]]
  delta_slope1 <- x[["delta_sl1"]]
  delta_slope2 <- x[["delta_sl2"]]
  
  lik <- pmap_dbl(
    my_tbl[, c("rt", "resp_recode", "pred_lr1", "pred_lr2")], 
    ~ dwiener(
      q = ..1, alpha = alpha, tau = tau, beta = beta, 
      delta = delta_ic + delta_slope1 * ..3 + delta_slope2 * ..4,
      resp = ..2
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


results_c_size <- optim(c(1, .1, .1, 1, .05), wiener_reg1_delta_log, my_tbl = tbl_cd[1:1000, ] %>% mutate(pred_lr = c_size_abs))


# now use these parameters to predict rt data and choices for the same stimuli used in the two-dimensional categorization task
l_tbl_up_and_down <- readRDS("data/l_tbl_up_and_down.RDS")
tbl_hotspots <- readRDS("data/hotspot-data.RDS")
tbl_transfer <- readRDS("data/transfer-rt-data.RDS")
tbl_test <- rbind(
  tbl_hotspots %>% dplyr::select(x1, x2, category, trial_id) %>% mutate(label = "old"),
  tbl_transfer  %>% dplyr::select(x1, x2, category, trial_id) %>% mutate(label = "new")
)

params_fin <- c(1.23, .5, .5)
lo <- c(0, .01, .0001, 1e-10)
hi <- c(10, .99, .9999, .999999999)

params <- pmap(list(params_fin[1:3], lo[1:3], hi[1:3]), upper_and_lower_bounds)

sims_hotspots <- pmap_dbl(
  tbl_test[, c("x1", "x2")], 
  ~ sum(pmap_dbl(
    tbl_hotspots[, c("x1", "x2")], 
    f_similarity, 
    c(params_fin[[2]], 1 - params_fin[[2]]), 15, tibble(x1 = .x, x2 = .y), 1
  ))
)

sims_strat <- pmap_dbl(
  tbl_test[, c("x1", "x2")], 
  ~ sum(pmap_dbl(
    l_tbl_up_and_down[[2]][, c("x1", "x2")], 
    f_similarity, 
    c(params_fin[[2]], 1 - params_fin[[2]]), 15, tibble(x1 = .x, x2 = .y), 1
  ))
)
cor(sims_strat, sims_hotspots)

p_gamma <- .2
p_thx_response <- .2
sims_strat_z <- scale(sims_strat^p_gamma)[, 1] #/ max(sims_strat^p_gamma)
sims_hotspots_z <- scale(sims_hotspots^p_gamma)[, 1] # / max(sims_hotspots^p_gamma)
hist(sims_strat_z)
hist(sims_hotspots_z)

cor(sims_strat_z, sims_hotspots_z)

plot(sort(sims_strat)/(sort(sims_strat) + p_thx_response))
plot(sort(sims_hotspots)/(sort(sims_hotspots) + p_thx_response))


tbl_test$sim_strat_z <- sims_strat_z
tbl_test$sim_hotspots_z <- sims_hotspots_z
tbl_test$pred_diff <- tbl_test$sim_strat_z - tbl_test$sim_hotspots_z
tbl_test$pred_diff_abs <- abs(tbl_test$pred_diff)
tbl_test$disc_power <- rank(desc(tbl_test$pred_diff_abs))

tbl_test_disc <- tbl_test %>% filter(disc_power <= 10)

ggplot(tbl_test_disc, aes(x1, x2)) +
  geom_point(aes(size = pred_diff, color = pred_diff)) + 
  ggrepel::geom_label_repel(aes(label = round(pred_diff, 2))) +
  theme_bw() +
  scale_x_continuous(expand = c(.03, 0)) +
  scale_y_continuous(expand = c(.03, 0)) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme(strip.background = element_rect(fill = "white")) +
  scale_size_continuous(name = "Pred. Difference", range = c(1, 5), guide = "none") +
  scale_color_gradient2(midpoint = 0, low = "skyblue2", high = "tomato4", guide = "none") +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 7))

n_reps <- 20

#saveRDS(tbl_transfer, file = "data/transfer-rt-data.RDS")

tbl_rt_strat <- map_df(
  tbl_test$sim_strat_z, 
  ~ rwiener(
    n = n_reps, alpha = results_c_size$par[[1]], tau = results_c_size$par[[2]], 
    beta = results_c_size$par[[3]], delta = results_c_size$par[[4]] + 2 *.x#results_c_size$par[[5]] * .x
  )
) %>% as_tibble() %>% 
  rename(rt = q) %>%
  mutate(
    model = "Strat. Sampling",
    resp_recode = resp,
    resp = fct_relabel(resp, ~ c("old", "new")),
    sim_strat_z = rep(tbl_test$sim_strat_z, each = n_reps),
    sim_hotspot_z = rep(tbl_test$sim_hotspots_z, each = n_reps)
  )

tbl_rt_hotspots <- map_df(
  tbl_test$sim_hotspots_z, 
  ~ rwiener(
    n = n_reps, alpha = results_c_size$par[[1]], tau = results_c_size$par[[2]], 
    beta = results_c_size$par[[3]], delta = results_c_size$par[[4]] + 2 *.x#results_c_size$par[[5]] * .x
  )
) %>% as_tibble() %>% 
  rename(rt = q) %>%
  mutate(
    model = "Hotspot",
    resp_recode = resp,
    resp = fct_relabel(resp, ~ c("old", "new")),
    sim_strat_z = rep(tbl_test$sim_strat_z, each = n_reps),
    sim_hotspot_z = rep(tbl_test$sim_hotspots_z, each = n_reps)
  )

tbl_test$item_id <- 1:nrow(tbl_test)
tbl_rt_hotspots$item_id <- rep(tbl_test$item_id, each = n_reps)
tbl_rt_strat$item_id <- rep(tbl_test$item_id, each = n_reps)
tbl_rt_both <- rbind(tbl_rt_strat, tbl_rt_hotspots)

tbl_acc <- grouped_agg(tbl_rt_both %>% mutate(p_old = resp == "old"), c(model, resp), c(rt)) %>%
  mutate(prop_resp = n / sum(n))

ggplot(tbl_rt_both, aes(rt)) +
  geom_histogram(fill = "skyblue2", color = "black") +
  geom_label(data = tbl_acc, aes(x = 1.5, y = n/10, label = str_c("Prop. = ", round(prop_resp, 3)))) +
  facet_grid(model ~ resp) + 
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "RT (s)", y = "Nr. Responses") +
  theme(strip.background = element_rect(fill = "white"))



tbl_test$label <- factor(tbl_test$label, labels = c(1, 0))
levels(tbl_rt_hotspots$resp) <- c(0, 1)
levels(tbl_rt_strat$resp) <- c(0, 1)

tbl_hotspots_agg <- tbl_rt_hotspots %>% 
  group_by(item_id) %>%
  summarize(p_new = mean(as.numeric(as.character(resp))))
tbl_strat_agg <- tbl_rt_strat %>%
  group_by(item_id) %>%
  summarize(p_new = mean(as.numeric(as.character(resp))))

hist(tbl_strat_agg$p_new)
hist(tbl_hotspots_agg$p_new)

plot(tbl_strat_agg$p_new, tbl_hotspots_agg$p_new)

tbl_test$p_new_hotspots <- tbl_hotspots_agg$p_new
tbl_test$p_new_strat <- tbl_strat_agg$p_new

# okish for recovery
tbl_test %>% group_by(label) %>%
  summarize(
    p_new_hotspots = mean(p_new_hotspots),
    p_new_strat = mean(p_new_strat)
  )

# map2 over sim_strat_z and sim_hotspot_z at the same time
# assume a larger coefficient on sim_hotspot_z than on sim_strat_z

# recover by model just assuming sim_hotspot_z as a covariate
# adapt wiener_ll function to handle to covariates
tbl_ll_diff <- tibble(
  n_reps = numeric(), 
  it = numeric(), 
  ll_diff = numeric(),
  ll_bivar = numeric(),
  ll_hotspot = numeric(),
  ll_strat = numeric()
)

for (nr in c(1, 2)) {
  for (i in 1:10) {
    tbl_rt_strat <- map2_df(
      tbl_test$sim_hotspots_z, tbl_test$sim_strat_z, 
      ~ rwiener(
        n = nr, alpha = results_c_size$par[[1]], 
        tau = results_c_size$par[[2]], beta = results_c_size$par[[3]], 
        delta = results_c_size$par[[4]] + 3 * .x  + 1 * .y# results_c_size$par[[4]]
      )
    ) %>% as_tibble() %>% 
      rename(rt = q) %>%
      mutate(
        model = "Strat. Sampling",
        resp_recode = resp,
        resp = fct_relabel(resp, ~ c("old", "new")),
        sim_strat_z = rep(tbl_test$sim_strat_z, each = nr),
        sim_hotspot_z = rep(tbl_test$sim_hotspots_z, each = nr)
      )
    
    results_recover_bivar_by_bivar <- optim(
      c(1, .1, .1, 1, 1, 1), 
      wiener_reg2_delta_log, 
      my_tbl = tbl_rt_strat %>% mutate(pred_lr1 = sim_hotspot_z, pred_lr2 = sim_strat_z)
    )
    results_recover_bivar_by_hotspot <- optim(
      c(1, .1, .1, 1, 1), 
      wiener_reg1_delta_log, 
      my_tbl = tbl_rt_strat %>% mutate(pred_lr = sim_hotspot_z)
    )
    results_recover_bivar_by_strat <- optim(
      c(1, .1, .1, 1, 1), 
      wiener_reg1_delta_log, 
      my_tbl = tbl_rt_strat %>% mutate(pred_lr = sim_strat_z)
    )
    tbl_ll_diff <- rbind(
      tbl_ll_diff,
      tibble(
        n_reps = nr,
        it = i,
        ll_bivar = results_recover_bivar_by_bivar$value,
        ll_hotspot = results_recover_bivar_by_hotspot$value,
        ll_strat = results_recover_bivar_by_strat$value
      )
    )
    saveRDS(tbl_ll_diff, file = "data/wiener-strat-sampling-recovery.RDS")
  }
}


tbl_ll_diff <- readRDS("data/wiener-strat-sampling-recovery.RDS")

tbl_ll_diff_long <- tbl_ll_diff %>% 
  mutate(
    n_data = n_reps * 330,
    bic_bivar = ll_bivar + 6*log(n_data), 
    bic_hotspot = ll_hotspot + 5*log(n_data), 
    bic_strat = ll_strat + 5*log(n_data), 
    aic_bivar = ll_bivar + 12,
    aic_hotspot = ll_hotspot + 10,
    aic_strat = ll_strat + 10,
    is_recovered_ll = ll_bivar == pmin(ll_bivar, ll_hotspot, ll_strat),
    is_recovered_bic = bic_bivar == pmin(bic_bivar, bic_hotspot, bic_strat),
    is_recovered_aic = aic_bivar == pmin(aic_bivar, aic_hotspot, aic_strat),
  ) %>%
  #pivot_longer(c(bic_bivar, bic_hotspot, bic_strat)) %>%
  pivot_longer(c(aic_bivar, aic_hotspot, aic_strat)) %>%
  #pivot_longer(c(ll_bivar, ll_hotspot, ll_strat)) %>%
  group_by(n_reps, it) %>% mutate(rwn = row_number(it)) %>% ungroup()

tbl_recovered_agg <- grouped_agg(
  tbl_ll_diff_long %>% filter(rwn == 1), 
  c(name, n_reps), 
  c(is_recovered_ll, is_recovered_bic, is_recovered_aic)
)

ggplot(tbl_ll_diff_long, aes(value, group = name)) +
  geom_histogram() +
  geom_label(
    data = tbl_recovered_agg, aes(
      x = mean(tbl_ll_diff_long$value), 
      y = tbl_recovered_agg$n[1]/4, 
      label = str_c("p (recov.) = ", round(mean_is_recovered_aic, 2))
    )
  ) +
  facet_grid(name ~ n_reps)


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




