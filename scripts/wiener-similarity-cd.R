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

tbl_cd <- tbl_cd[sample(1:nrow(tbl_cd), nrow(tbl_cd), replace = FALSE) < (nrow(tbl_cd) / 8), ]

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

set.seed(0)
dat <- rwiener(n=100, alpha=2, tau=.3, beta=.5, delta=.5)
optim1 <- optim(c(1, .1, .1, 1), wiener_deviance, dat=dat, method="Nelder-Mead")

dat <- rwiener(n=10000, alpha=2, tau=.3, beta=.5, delta=.5)
optim2 <- optim(c(1, .1, .1, 1), wiener_deviance, dat=dat, method="Nelder-Mead")


dwiener(.75, 2, .3, .5, .5)
plot(dwiener(seq(0, 5, by = .01), 2, .3, .5, .5, resp = rep("upper", 501)))



wiener_reg_delta_log <- function(x, my_tbl) {
  alpha <- x[[1]]
  tau <- x[[2]]
  beta <- x[[3]]
  delta_ic <- x[[4]]
  delta_slope <- x[[5]]
  
  lik <- pmap_dbl(
    my_tbl[, c("rt", "resp_recode", "setSize")], ~ dwiener(
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

tbl_cd$resp_recode <- map_chr(tbl_cd$resp + 1, ~ c("upper", "lower")[.x])

optim(c(1, .1, .1, 1, .05), wiener_reg_delta_log, my_tbl = tbl_cd[1:1000, ])








