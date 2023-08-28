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
  alpha = 1.3, beta = .55, delta_ic = .9, delta_sl1 = 3, delta_sl2 = 1, tau = .25
)
x <- params_gen
w <- params_gen[["w"]]
sens <- params_gen[["sens"]] 
gamma <- params_gen[["gamma"]] 
sim_hs <- l_sims$sims_hotspots_z
sim_strat <- l_sims$sims_strat_z
alpha <- params_gen[["alpha"]]
beta <- params_gen[["beta"]]
tau <- params_gen[["tau"]]
delta_ic <- params_gen[["delta_ic"]]
delta_sl1 <- params_gen[["delta_sl1"]]
delta_sl2 <- params_gen[["delta_sl2"]]
n_reps <- 1


tbl_gen <- gen_recognition_2_slopes(params_gen, tbl_hotspots, tbl_strat, tbl_test)

# testing
ll_recognition_2_slopes(
  params_gen, tbl_hotspots, tbl_strat, 
  tbl_gen[, c("x1", "x2", "resp_recode", "rt")]
)

ll_recognition_1_slope(
  params_gen, tbl_hotspots, tbl_strat, 
  tbl_gen[, c("x1", "x2", "resp_recode", "rt")], "sims_hotspots_z"
)

ll_recognition_1_slope(
  params_gen, tbl_hotspots, tbl_strat, 
  tbl_gen[, c("x1", "x2", "resp_recode", "rt")], "sims_strat_z"
)

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
    method = "Nelder-Mead", control = list(maxit = 2000)
  )
  results_recover_hotspot <- optim(
    params_init_1, ll_recognition_1_slope, 
    tbl_hotspots = tbl_hotspots, 
    tbl_strat = tbl_strat, 
    tbl_test = tbl_gen[, c("x1", "x2", "resp_recode", "rt")],
    sim_which = "sims_hotspots_z",
    method = "Nelder-Mead", control = list(maxit = 2000)
  )
  results_recover_strat <- optim(
    params_init_1, ll_recognition_1_slope, 
    tbl_hotspots = tbl_hotspots, 
    tbl_strat = tbl_strat, 
    tbl_test = tbl_gen[, c("x1", "x2", "resp_recode", "rt")],
    sim_which = "sims_strat_z",
    method = "Nelder-Mead", control = list(maxit = 2000)
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
  1:10, recovery_study,
  tbl_hotspots = tbl_hotspots,
  tbl_strat = tbl_strat,
  tbl_test = tbl_test,
  params_init_1 = params_init_1,
  params_init_2 = params_init_2,
  params_gen = params_gen,
  .progress = TRUE, .options = furrr_options(seed = TRUE)
)  





