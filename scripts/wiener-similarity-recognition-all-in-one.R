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
conflicts_prefer(dplyr::filter)


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
  w = .5, sens = 2, gamma = .2, 
  alpha = 1.3, beta = .55, tau = .25, 
  delta_ic = .9, delta_sl1 = 3, 
  # the following two parameters are optional
  delta_sl2 = 1,
  sens2 = 4
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

tbl_gen_1_sens_1_slope <- gen_recognition(params_gen_1_sens_1_slope, tbl_train_not_important, tbl_important, tbl_recognition)[[1]]
tbl_gen_2_sens_1_slope <- gen_recognition(params_gen_2_sens_1_slope, tbl_train_not_important, tbl_important, tbl_recognition)[[1]]
tbl_gen_1_sens_2_slopes <- gen_recognition(params_gen_1_sens_2_slopes, tbl_train_not_important, tbl_important, tbl_recognition)[[1]]
tbl_gen_2_sens_2_slopes <- gen_recognition(params_gen_2_sens_2_slopes, tbl_train_not_important, tbl_important, tbl_recognition)[[1]]

# plausibility checks
tbl_gen_1_sens_1_slope %>%
  group_by(label) %>%
  count(resp)

tbl_gen_1_sens_2_slopes %>%
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

lo <- c(.01, 0, .01, 4, .9, .2, 5, 5, 5, 0)
hi <- c(.99, 40, .99, .1, .1, .4, .01, .01, .01, 40)
params_init_2_2 <- pmap(list(start_vals_2, lo, hi), upper_and_lower_bounds)
# notation _1_2 means one sensitivity and two variables regressed on drift rate



# testing

tbl_gen_2_sens_1_slopes <- gen_recognition(params_gen_2_sens_1_slope, tbl_train_not_important, tbl_important, tbl_recognition)[[1]]


ll_recognition_1_1(
  pmap(list(params_gen_1_sens_1_slope, as.list(lo[1:8]), as.list(hi[1:8])), upper_and_lower_bounds),
  tbl_train_not_important,
  tbl_important,
  tbl_gen_2_sens_1_slopes[, c("x1", "x2", "resp_recode", "rt")],
  lo[1:8], hi[1:8]
)


ll_recognition_1_2(
  pmap(list(params_gen_1_sens_2_slopes, as.list(lo[1:9]), as.list(hi[1:9])), upper_and_lower_bounds),
  tbl_train_not_important,
  tbl_important,
  tbl_gen_2_sens_1_slopes[, c("x1", "x2", "resp_recode", "rt")],
  lo[1:9], hi[1:9]
)


ll_recognition_2_1(
  pmap(list(params_gen_2_sens_1_slope, as.list(lo[c(1:8, 10)]), as.list(hi[c(1:8, 10)])), upper_and_lower_bounds),
  tbl_train_not_important,
  tbl_important,
  tbl_gen_2_sens_1_slopes[, c("x1", "x2", "resp_recode", "rt")],
  lo[c(1:8, 10)], hi[c(1:8, 10)]
)


ll_recognition_2_2(
  pmap(list(params_gen_2_sens_2_slopes, as.list(lo), as.list(hi)), upper_and_lower_bounds),
  tbl_train_not_important,
  tbl_important,
  tbl_gen_2_sens_1_slopes[, c("x1", "x2", "resp_recode", "rt")],
  lo, hi
)




recovery_study <- function(
    i, n_reps, tbl_train, tbl_important, tbl_test, 
    params_init_2_2, params_gen,
    string_save
) {
  params_init_1_1 <- params_init_2_2[1:8]
  params_init_1_2 <- params_init_2_2[1:9]
  params_init_2_1 <- params_init_2_2[c(1:8, 10)]
  
  tbl_recognition_rep <- repeat_tibble(tbl_recognition, n_reps)
  tbl_gen <- gen_recognition(params_gen, tbl_train, tbl_important, tbl_recognition_rep)
  results_recover_1_1 <- optim(
    params_init_1_1, ll_recognition_1_1, 
    tbl_train = tbl_train, 
    tbl_strat = tbl_important, 
    tbl_test = tbl_gen[[1]][, c("x1", "x2", "resp_recode", "rt")],
    lo = lo[1:8],
    hi = hi[1:8],
    method = "Nelder-Mead", control = list(maxit = 1000)
  )
  results_recover_1_2 <- optim(
    params_init_1_2, ll_recognition_1_2, 
    tbl_train = tbl_train, 
    tbl_strat = tbl_important, 
    tbl_test = tbl_gen[[1]][, c("x1", "x2", "resp_recode", "rt")],
    lo = lo[1:9],
    hi = hi[1:9],
    method = "Nelder-Mead", control = list(maxit = 1000)
  )
  results_recover_2_1 <- optim(
    params_init_2_1, ll_recognition_2_1, 
    tbl_train = tbl_train, 
    tbl_strat = tbl_important, 
    tbl_test = tbl_gen[[1]][, c("x1", "x2", "resp_recode", "rt")],
    lo = lo[c(1:8, 10)],
    hi = hi[c(1:8, 10)],
    method = "Nelder-Mead", control = list(maxit = 1000)
  )
  results_recover_2_2 <- optim(
    params_init_2_2, ll_recognition_2_2, 
    tbl_train = tbl_train, 
    tbl_strat = tbl_important, 
    tbl_test = tbl_gen[[1]][, c("x1", "x2", "resp_recode", "rt")],
    lo = lo,
    hi = hi,
    method = "Nelder-Mead", control = list(maxit = 1000)
  )
  l_out <- list(
    it = i,
    '1_1' = results_recover_1_1,
    '1_2' = results_recover_1_2,
    '2_1' = results_recover_2_1,
    '2_2' = results_recover_2_2
  )
  saveRDS(l_out, file = str_c(
    "data/representation-and-decision/recovery-representation-and-decision-", string_save, "-", i, ".RDS"
  ))
  
  return(l_out)
}



# one sensitivity and, but two similarities affecting drift rate differently
n_datasets <- 4*14

# 1 1
t_start <- Sys.time()
plan(multisession, workers = future::availableCores() - 2)
l_recovery_results_1_1 <- future_map2(
  1:n_datasets, rep(1, n_datasets) , recovery_study,
  tbl_train = tbl_train_not_important,
  tbl_important = tbl_important,
  tbl_test = tbl_recognition,
  params_init_2_2 = params_init_2_2,
  params_gen = params_gen_1_sens_1_slope,
  string_save = "gen-1-1",
  .progress = TRUE, .options = furrr_options(seed = TRUE)
)
future::plan("default")
beepr::beep()
t_end <- Sys.time()
round(t_end - t_start, 1)

# 1 2
t_start <- Sys.time()
plan(multisession, workers = future::availableCores() - 2)
l_recovery_results_1_2 <- future_map2(
  1:n_datasets, rep(1, n_datasets) , recovery_study,
  tbl_train = tbl_train_not_important,
  tbl_important = tbl_important,
  tbl_test = tbl_recognition,
  params_init_2_2 = params_init_2_2,
  params_gen = params_gen_1_sens_2_slopes,
  string_save = "gen-1-2",
  .progress = TRUE, .options = furrr_options(seed = TRUE)
)
future::plan("default")
beepr::beep()
t_end <- Sys.time()
round(t_end - t_start, 1)

# 2 1
t_start <- Sys.time()
plan(multisession, workers = future::availableCores() - 2)
l_recovery_results_2_1 <- future_map2(
  1:n_datasets, rep(1, n_datasets) , recovery_study,
  tbl_train = tbl_train_not_important,
  tbl_important = tbl_important,
  tbl_test = tbl_recognition,
  params_init_2_2 = params_init_2_2,
  params_gen = params_gen_2_sens_1_slope,
  string_save = "gen-2-1",
  .progress = TRUE, .options = furrr_options(seed = TRUE)
)
future::plan("default")
beepr::beep()
t_end <- Sys.time()
round(t_end - t_start, 1)

# 2 2
t_start <- Sys.time()
plan(multisession, workers = future::availableCores() - 2)
l_recovery_results_2_2 <- future_map2(
  1:n_datasets, rep(1, n_datasets) , recovery_study,
  tbl_train = tbl_train_not_important,
  tbl_important = tbl_important,
  tbl_test = tbl_recognition,
  params_init_2_2 = params_init_2_2,
  params_gen = params_gen_2_sens_2_slopes,
  string_save = "gen-2-2",
  .progress = TRUE, .options = furrr_options(seed = TRUE)
)
future::plan("default")
beepr::beep()
t_end <- Sys.time()
round(t_end - t_start, 1)





file_paths <- str_c("data/representation-and-decision/", dir("data/representation-and-decision/"))
l_recovery_results <- map(file_paths, readRDS)

ll_1_1 <- map_dbl(map(l_recovery_results_1_1, "1_1"), "value")
ll_1_2 <- map_dbl(map(l_recovery_results_1_1, "1_2"), "value")
ll_2_1 <- map_dbl(map(l_recovery_results_1_1, "2_1"), "value")
ll_2_2 <- map_dbl(map(l_recovery_results_1_1, "2_2"), "value")

tbl_lls <- as_tibble(cbind(ll_1_1, ll_1_2, ll_2_1, ll_2_2))
tbl_lls$it <- 1:nrow(tbl_lls)

tbl_lls <- tbl_lls %>% mutate(
  n_data = c(rep(122, 14), rep(244, 14)),
  bic_1_1 = ll_1_1 + 8*log(n_data), 
  bic_1_2 = ll_1_2 + 9*log(n_data), 
  bic_2_1 = ll_2_1 + 9*log(n_data), 
  bic_2_2 = ll_2_2 + 10*log(n_data), 
  aic_1_1 = ll_1_1 + 16, 
  aic_1_2 = ll_1_2 + 18, 
  aic_2_1 = ll_2_1 + 18, 
  aic_2_2 = ll_2_2 + 20, 
  bic_min = pmin(bic_1_1, bic_1_2, bic_2_1, bic_2_2),
  is_recovered_ll = ll_1_1 == pmin(ll_1_1, ll_1_2, ll_2_1, ll_2_2),
  is_recovered_bic = bic_1_1 == bic_min,
  is_recovered_aic = aic_1_1 == pmin(aic_1_1, aic_1_2, aic_2_1, aic_2_2)
) 



tbl_lls_long <- tbl_lls %>%
  select(it, n_data, bic_1_1, bic_1_2, bic_2_1, bic_2_2, is_recovered_bic, bic_min) %>%
  pivot_longer(c(bic_1_1, bic_1_2, bic_2_1, bic_2_2)) %>%
  #pivot_longer(c(aic_bivar, aic_hotspot, aic_strat)) %>%
  #pivot_longer(c(ll_bivar, ll_hotspot, ll_strat)) %>%
  mutate(name = factor(name))

levels(tbl_lls_long$name) <- c("1 Sim, 1 Drift", "1 Sim, 2 Drifts", "2 Sims, 1 Drift", "2 Sims, 2 Drifts")

tbl_recovered_agg <- grouped_agg(
  tbl_lls_long %>% filter(bic_min == value), 
  c(name, n_data), #  
  c(is_recovered_bic)
) %>% select(name, n_data, n) %>%
  mutate(model_in = "1 Sim, 1 Drift")

ggplot(tbl_recovered_agg, aes(model_in, name)) +
  geom_tile(aes(fill = n), show.legend = FALSE) +
  geom_label(aes(label = n)) +
  facet_wrap(~ n_data) + 
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Model In", y = "Model Out") +
  theme(
    strip.background = element_rect(fill = "white"),
    text = element_text(size = 16)
  )


ggplot(tbl_lls_long, aes(value, group = name)) +
  geom_histogram(fill = "skyblue2", color = "black") +
  geom_label(
    data = tbl_recovered_agg, aes(
      x = mean(tbl_lls_long$value), 
      y = tbl_recovered_agg$n[1]/6, 
      label = str_c("p (recov.) = ", round(mean_is_recovered_bic, 2))
    )
  ) + facet_wrap(~ name) + # facet_grid(name ~ n_reps) +
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


