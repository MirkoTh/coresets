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

l_tbl_up_and_down <- readRDS("data/l_tbl_up_and_down.RDS")
tbl_hotspots <- readRDS("data/hotspot-data.RDS")
tbl_transfer <- readRDS("data/transfer-rt-data.RDS")
tbl_test <- rbind(
  tbl_hotspots %>% dplyr::select(x1, x2, category, trial_id) %>% mutate(label = "old"),
  tbl_transfer  %>% dplyr::select(x1, x2, category, trial_id) %>% mutate(label = "new")
)
tbl_strat <- l_tbl_up_and_down[[2]]

params_gen <- list(
  w = .5, sens = 15, 
  gamma = .2, 
  alpha = 1.3, beta = .55, tau = .25, delta_ic = .9, delta_sl1 = 3, delta_sl2 = 1
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
n_reps <- 1


tbl_gen <- gen_recognition_2_slopes(params_gen, tbl_hotspots, tbl_strat, tbl_test)

start_vals_2 <- list(
  w = .5,
  sens = 10,
  gamma = .2,
  alpha = 1,
  beta = .5,
  tau = .3,
  delta_ic = 1,
  delta_sl1 = 2,
  delta_sl2 = 2
)

lo2 <- c(.01, 0, .01, 4, .9, .2, 5, 5, 5)
hi2 <- c(.99, 20, .99, .1, .1, .4, .01, .01, .01)
params <- pmap(list(start_vals_2, lo2, hi2), upper_and_lower_bounds)
params_init_2 <- params
params_init_1 <- params[1:8]
lo1 <- lo2[1:8]
hi1 <- hi2[1:8]


params_gen_sc <- pmap(list(params_gen, lo2, hi2), upper_and_lower_bounds)

# testing
ll_recognition_2_slopes(
  params_gen, tbl_hotspots, tbl_strat, 
  tbl_gen[, c("x1", "x2", "resp_recode", "rt")]
)

ll_recognition_1_slope(
  params_gen[1:8], tbl_hotspots, tbl_strat, 
  tbl_gen[, c("x1", "x2", "resp_recode", "rt")], "sims_hotspots_z"
)

ll_recognition_1_slope(
  params_gen[1:8], tbl_hotspots, tbl_strat, 
  tbl_gen[, c("x1", "x2", "resp_recode", "rt")], "sims_strat_z"
)



recovery_study <- function(
    i, tbl_hotspots, tbl_strat, tbl_test, 
    params_init_1, params_init_2, params_gen
) {
  tbl_gen <- gen_recognition_2_slopes(params_gen, tbl_hotspots, tbl_strat, tbl_test)
  results_recover_bivar <- optim(
    params_init_2, ll_recognition_2_slopes, 
    tbl_hotspots = tbl_hotspots, 
    tbl_strat = tbl_strat, 
    tbl_test = tbl_gen[, c("x1", "x2", "resp_recode", "rt")],
    method = "Nelder-Mead", control = list(maxit = 1000)
  )
  results_recover_hotspot <- optim(
    params_init_1, ll_recognition_1_slope, 
    tbl_hotspots = tbl_hotspots, 
    tbl_strat = tbl_strat, 
    tbl_test = tbl_gen[, c("x1", "x2", "resp_recode", "rt")],
    sim_which = "sims_hotspots_z",
    method = "Nelder-Mead", control = list(maxit = 1000)
  )
  results_recover_strat <- optim(
    params_init_1, ll_recognition_1_slope, 
    tbl_hotspots = tbl_hotspots, 
    tbl_strat = tbl_strat, 
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
  tbl_strat = tbl_strat,
  tbl_test = tbl_test,
  params_init_1 = params_init_1,
  params_init_2 = params_init_2,
  params_gen = params_gen,
  .progress = TRUE, .options = furrr_options(seed = TRUE)
)  
l_recovery_results[[1]]


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
      y = tbl_recovered_agg$n[1]/4, 
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


