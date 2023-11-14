# only recognition, not categorization

# we have
# stimuli shown during training -> hotspot model
# stimuli strategically sampled -> strat model

# parameters to generate the data
# 2: sensitivity, attentional weighting
# 1: gamma to scale 
# 6: alpha, beta, delta_ic, delta_sl1, delta_sl2, tau

# 9 parameter generating model
# 8/9 parameter models to fit the generated data

# assumptions
# one set of w and sensitivity parameters (c)
# i.e., generalization based on presented and imagined data works similarly

# setup
#1 generate data
## calculate similarities given bias and sens
## generate choices and rts with wiener likelihood given similarities

#2 fit data
## calculate similarities given bias and sens
## get ll of choices and rts with wiener likelihood given similarities
## either hotspot similarities or strat similarities (competitor)
## or hotspot & strat similarities (generating)


rm(list = ls())
set.seed(3453)

library(conflicted)
library(tidyverse)
library(grid)
library(gridExtra)
library(furrr)
library(docstring)
library(rutils)
library(cmdstanr)
library(RWiener)
conflicts_prefer(RWiener::rwiener)


path_load <- c("utils/utils.R", "scripts/stan-wiener.R")
walk(path_load, source)





# all data seen during training
tbl_train <- readRDS("data/train-test-third-dimension.RDS")$train
# data presented during transfer without category feedback
tbl_transfer <-  readRDS("data/train-test-third-dimension.RDS")$transfer
tbl_recognition <- rbind(tbl_train %>% mutate(label = "old"), tbl_transfer %>% mutate(label = "new"))

l_tbl_important_down_0 <- readRDS(file = "data/downsampling-ii-uniform-cat-0-third-dim.RDS")
l_tbl_important_down_1 <- readRDS(file = "data/downsampling-ii-uniform-cat-1-third-dim.RDS")
l_tbl_important <- map(
  1:length(l_tbl_important_down_0), 
  l_0 = l_tbl_important_down_0, 
  l_1 = l_tbl_important_down_1,
  extract_n_most_important_points
)
# exemplary set with average number of retained data points
tbl_important <- l_tbl_important[[4]]
# throw out important data points from separate training data set
tbl_train_not_important <- tbl_train %>% 
  left_join(tbl_important[, c("x1", "x2", "rank_importance")]) %>%
  filter(is.na(rank_importance)) %>%
  select(-rank_importance)

params_gen <- list(
  w = .5, sens = 15, gamma = .2, 
  alpha = 1.3, beta = .55, tau = .25, 
  delta_ic = .9, delta_sl1 = 3, 
  # the following two parameters are optional
  delta_sl2 = 1,
  sens2 = 20
)

x <- params_gen
w <- params_gen[["w"]]
sens <- params_gen[["sens"]] 
gamma <- params_gen[["gamma"]] 
alpha <- params_gen[["alpha"]]
beta <- params_gen[["beta"]]
tau <- params_gen[["tau"]]
delta_ic <- params_gen[["delta_ic"]]
delta_sl1 <- params_gen[["delta_sl1"]]
delta_sl2 <- params_gen[["delta_sl2"]]
sens2 <- params_gen[["sens2"]]
n_reps <- 1

params_gen_1_sens_1_slope <- params_gen[!names(params_gen) %in% c("delta_sl2", "sens2")]
params_gen_2_sens_1_slope <- params_gen[!names(params_gen) %in% c("delta_sl2")]
params_gen_1_sens_2_slopes <- params_gen[!names(params_gen) %in% c("sens2")]
params_gen_2_sens_2_slopes <- params_gen


tbl_gen_1_sens_1_slope <- gen_recognition(params_gen_1_sens_1_slope, tbl_train_not_important, tbl_important, tbl_recognition)
tbl_gen_2_sens_1_slope <- gen_recognition(params_gen_2_sens_1_slope, tbl_train_not_important, tbl_important, tbl_recognition)
tbl_gen_1_sens_2_slopes <- gen_recognition(params_gen_1_sens_2_slopes, tbl_train_not_important, tbl_important, tbl_recognition)
tbl_gen_2_sens_2_slopes <- gen_recognition(params_gen_2_sens_2_slopes, tbl_train_not_important, tbl_important, tbl_recognition)

tbl_gen_1_sens_1_slope %>%
  group_by(label) %>%
  count(resp)

tbl_gen_2_sens_1_slope %>%
  group_by(label) %>%
  count(resp)

ggplot(tbl_gen_2_sens_1_slope, aes(rt)) +
  geom_histogram() +
  facet_grid(label ~ resp)

start_vals_2 <- list(
  w = .5,
  sens = 10,
  gamma = .2,
  alpha = 1,
  beta = .5,
  tau = .3,
  delta_ic = 1,
  delta_sl1 = 2,
  delta_sl2 = 2,
  sens2 = 15
)

lo_2_2 <- c(.01, 0, .01, 4, .9, .2, 5, 5, 5, 0)
hi_2_2 <- c(.99, 40, .99, .1, .1, .4, .01, .01, .01, 40)
params <- pmap(list(start_vals_2, lo_2_2, hi_2_2), upper_and_lower_bounds)
# notation _1_2 means one sensitivity and two variables regressed on drift rate
params_init_2_2 <- params
params_init_1_1 <- params[1:8]
params_init_2_1 <- params[c(1:8, 10)]
params_init_1_2 <- params[1:9]
lo_1_1 <- lo_2_2[1:8]
hi_1_1 <- hi_2_2[1:8]
lo_2_1 <- lo_2_2[c(1:8, 10)]
hi_2_1 <- hi_2_2[c(1:8, 10)]
lo_1_2 <- lo_2_2[1:9]
hi_1_2 <- hi_1_2[1:9]


params_gen_sc <- pmap(list(params_gen, lo2, hi2), upper_and_lower_bounds)

# testing
ll_recognition(
  params_gen_2_sens_1_slope,
  tbl_train_not_important,
  tbl_important,
  tbl_gen_2_sens_1_slope[, c("x1", "x2", "resp_recode", "rt")],
  lo_2_1, hi_2_1
)


ll_recognition_2_slopes(
  params_gen_2_sens_2_slopes, tbl_train, tbl_important, 
  tbl_gen_2[, c("x1", "x2", "resp_recode", "rt")]
)

ll_recognition_1_slope(
  params_gen[1:8], tbl_train, tbl_important, 
  tbl_gen_2[, c("x1", "x2", "resp_recode", "rt")], "sims_hotspots_z"
)

ll_recognition_1_slope(
  params_gen[1:8], tbl_train, tbl_important, 
  tbl_gen_2[, c("x1", "x2", "resp_recode", "rt")], "sims_strat_z"
)



recovery_study <- function(
    i, tbl_hotspots, tbl_important, tbl_test, 
    params_init_1, params_init_2, params_gen
) {
  tbl_gen <- gen_recognition(params_gen, tbl_train, tbl_important, tbl_recognition)
  results_recover_bivar <- optim(
    params_init_2, ll_recognition_2_slopes, 
    tbl_hotspots = tbl_hotspots, 
    tbl_strat = tbl_important, 
    tbl_test = tbl_gen[, c("x1", "x2", "resp_recode", "rt")],
    method = "Nelder-Mead", control = list(maxit = 1000)
  )
  results_recover_hotspot <- optim(
    params_init_1, ll_recognition_1_slope, 
    tbl_hotspots = tbl_hotspots, 
    tbl_strat = tbl_important, 
    tbl_test = tbl_gen[, c("x1", "x2", "resp_recode", "rt")],
    sim_which = "sims_hotspots_z",
    method = "Nelder-Mead", control = list(maxit = 1000)
  )
  results_recover_strat <- optim(
    params_init_1, ll_recognition_1_slope, 
    tbl_hotspots = tbl_hotspots, 
    tbl_strat = tbl_important, 
    tbl_test = tbl_gen[, c("x1", "x2", "resp_recode", "rt")],
    sim_which = "sims_strat_z",
    method = "Nelder-Mead", control = list(maxit = 1000)
  )
  l_out <- list(
    it = i,
    bivar = results_recover_bivar,
    hotspot = results_recover_hotspot,
    strat = results_recover_strat
  )
  saveRDS(l_out, file = str_c(
    "data/representation-and-decision/recovery-representation-and-decision-", i, ".RDS"
  ))
  
  return(l_out)
}



plan(multisession, workers = future::availableCores() - 2)

l_recovery_results <- future_map(
  1:60, recovery_study,
  tbl_hotspots = tbl_hotspots,
  tbl_strat = tbl_important,
  tbl_test = tbl_test,
  params_init_1 = params_init_1,
  params_init_2 = params_init_2,
  params_gen = params_gen,
  .progress = TRUE, .options = furrr_options(seed = TRUE)
)  

file_paths <- str_c("data/representation-and-decision/", dir("data/representation-and-decision/"))
l_recovery_results <- map(file_paths, readRDS)

ll_bivar <- map_dbl(map(l_recovery_results, "bivar"), "value")
ll_hotspot <- map_dbl(map(l_recovery_results, "hotspot"), "value")
ll_strat <- map_dbl(map(l_recovery_results, "strat"), "value")

tbl_lls <- as_tibble(cbind(ll_bivar, ll_hotspot, ll_strat))
tbl_lls$it <- 1:nrow(tbl_lls)

tbl_lls_long <- tbl_lls %>% mutate(
  n_reps = 1,
  n_data = n_reps * 330,
  bic_bivar = ll_bivar + 9*log(n_data), 
  bic_hotspot = ll_hotspot + 8*log(n_data), 
  bic_strat = ll_strat + 8*log(n_data), 
  aic_bivar = ll_bivar + 18,
  aic_hotspot = ll_hotspot + 16,
  aic_strat = ll_strat + 16,
  is_recovered_ll = ll_bivar == pmin(ll_bivar, ll_hotspot, ll_strat),
  is_recovered_bic = bic_bivar == pmin(bic_bivar, bic_hotspot, bic_strat),
  is_recovered_aic = aic_bivar == pmin(aic_bivar, aic_hotspot, aic_strat), 
) %>%
  pivot_longer(c(bic_bivar, bic_hotspot, bic_strat)) %>%
  #pivot_longer(c(aic_bivar, aic_hotspot, aic_strat)) %>%
  #pivot_longer(c(ll_bivar, ll_hotspot, ll_strat)) %>%
  group_by(n_reps, it) %>% mutate(rwn = row_number(it)) %>% ungroup() %>%
  mutate(name = factor(name))

levels(tbl_lls_long$name) <- c("Bivariate", "Hotspot", "Strategic")

tbl_recovered_agg <- grouped_agg(
  tbl_lls_long %>% filter(rwn == 1), 
  c(name, n_reps), 
  c(is_recovered_ll, is_recovered_bic, is_recovered_aic)
)

ggplot(tbl_lls_long, aes(value, group = name)) +
  geom_histogram(fill = "skyblue2", color = "black") +
  geom_label(
    data = tbl_recovered_agg, aes(
      x = mean(tbl_lls_long$value), 
      y = tbl_recovered_agg$n[1]/6, 
      label = str_c("p (recov.) = ", round(mean_is_recovered_bic, 2))
    )
  ) + facet_grid(name ~ n_reps) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "BIC", y = "Nr. Observations") + 
  theme(strip.background = element_rect(fill = "white"))


l_pars_bivar_sc <- map(map(l_recovery_results, "bivar"), "par")

l_pars <- map(l_pars_bivar_sc, ~ pmap(list(.x, lo2, hi2), upper_and_lower_bounds_revert))
l_pars <- map(l_pars, unlist)
tbl_pars <- as_tibble(as.data.frame(reduce(l_pars, rbind)))
tbl_pars$it <- 1:nrow(tbl_pars)

tbl_pars_gen <- as_tibble(params_gen) %>% mutate(it = 1) %>%
  pivot_longer(-it)
tbl_pars_start <- as_tibble(start_vals_2) %>% mutate(it = 1) %>%
  pivot_longer(-it)

tbl_pars %>% pivot_longer(-it) %>%
  ggplot(aes(value)) +
  geom_vline(data = tbl_pars_start, aes(xintercept = value), color = "red", size = 1.5) +
  geom_vline(data = tbl_pars_gen, aes(xintercept = value), color = "forestgreen", size = 1.5) +
  geom_histogram(fill = "skyblue2", color = "black") +
  facet_wrap(~ name, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Parameter Value", y = "Nr. Fits") + 
  theme(strip.background = element_rect(fill = "white"))


