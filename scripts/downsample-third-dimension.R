# Load packages and utils -------------------------------------------------
devtools::install_github("r-lib/conflicted")

rm(list = ls())
set.seed(43995)

library(conflicted)
library(docstring)
library(tidyverse)
conflicts_prefer(dplyr::filter())
library(grid)
library(gridExtra)
library(furrr)
library(rutils)
library(e1071)


path_load <- c("utils/utils.R")
walk(path_load, source)

is_fitting <- TRUE


# we want seven values of orientation and length
# people cannot discriminate > 7 values without training (e.g., Miller 1956)

# do not show same stimuli during training and transfer of categorization task

n_feat <- 2
d_measure <- 1
l_x1 <- list(train = seq(1, 7, by = 1), transfer = seq(.5, 7.5, by = 1))
l_x2 <- list(train = seq(.5, 7.5, by = 1), transfer = seq(1, 7, by = 1))

train_and_transfer_sets <- function(x1, x2) {
  tbl_x_ii <- crossing(x1, x2)
  tbl_x_ii$category <- factor(as.numeric(tbl_x_ii$x1 < tbl_x_ii$x2))
  tbl_x_ii <- tbl_x_ii[!(tbl_x_ii$x1 == tbl_x_ii$x2), ]
  tbl_x_ii$cat_structure <- "information-integration"
  tbl_x_ii$trial_id <- seq(1, nrow(tbl_x_ii), by = 1)
  # randomize trial order (only relevant for forgetful gcm)
  tbl_x_ii <- tbl_x_ii[
    sample(1:nrow(tbl_x_ii), nrow(tbl_x_ii), replace = FALSE), 
  ]
  tbl_x_ii$trial_id <- seq(1, nrow(tbl_x_ii), by = 1)
  return(tbl_x_ii)
}

# generate cl train and test sets
# simulate responses using assumed difficulty function
l_tbl_x_ii <- map2(l_x1, l_x2, train_and_transfer_sets) %>%
  map(., simulate_responses)
plot_grid(l_tbl_x_ii$train) + geom_abline()


if (is_fitting) {
  saveRDS(l_tbl_x_ii, file = "data/train-test-third-dimension.RDS")
} else if (!is_fitting) {
  l_tbl_x_ii <- readRDS(file = "data/train-test-third-dimension.RDS")
}

n_bins <- 5
l_tbl_x_ii$train %>% 
  mutate(x1_cut = cut(x1, n_bins), x2_cut = cut(x2, n_bins)) %>%
  group_by(x1_cut, x2_cut) %>%
  summarize(prop_correct = mean(p_correct)) %>%
  ungroup() %>%
  ggplot(aes(x1_cut, x2_cut)) +
  geom_tile(aes(fill = prop_correct)) +
  geom_text(aes(label = round(prop_correct, 2))) +
  scale_fill_gradient(name = "Accuracy", low = "white", high = "forestgreen", limits = c(0, 1)) + 
  geom_abline() + ggtitle("Simulated Response Accuracy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# fit plain-vanilla gcm on simulated cl training choices ------------------


params <- c(c = 1, w = .5, bias = .5, delta = .5)

lo <- c(0, .01, .0001, 1e-10)
hi <- c(10, .99, .9999, .999999999)
params_init <- pmap(list(params[1:3], lo[1:3], hi[1:3]), upper_and_lower_bounds)


# takes approx 2 min to run on laptop with 300 train and transfer stimuli
# takes approx 15 seconds to run on laptop with 105 train and transfer stimuli
t_start <- Sys.time()
results_pre <- optim(
  params_init,
  gcm_likelihood_no_forgetting,
  tbl_transfer = l_tbl_x_ii$train,
  tbl_x = l_tbl_x_ii$train, 
  n_feat = 2,
  d_measure = 1,
  lo = lo[1:3],
  hi = hi[1:3]
)
t_end <- Sys.time()
round(t_end - t_start, 1)

params_fin <- list()
params_fin[["not_tf"]] <- pmap_dbl(list(results_pre$par, lo[1:3], hi[1:3]), upper_and_lower_bounds_revert)
params_fin[["tf"]] <- results_pre$par



# fit forgetful gcm on simulated cl training choices ----------------------



params_init_forget <- pmap(list(params, lo, hi), upper_and_lower_bounds)


# takes approx 2 min to run on laptop with 300 train and transfer stimuli
# takes approx 15 seconds to run on laptop with 105 train and transfer stimuli
t_start <- Sys.time()
results_forget <- optim(
  params_init_forget,
  gcm_likelihood_forgetting,
  tbl_transfer = l_tbl_x_ii$train,
  tbl_x = l_tbl_x_ii$train, 
  n_feat = 2,
  d_measure = 1,
  lo = lo,
  hi = hi
)
t_end <- Sys.time()
round(t_end - t_start, 1)

params_fin_forget <- list()
params_fin_forget[["not_tf"]] <- pmap_dbl(list(results_forget$par, lo, hi), upper_and_lower_bounds_revert)
params_fin_forget[["tf"]] <- results_forget$par



# throw out less important points up to K ---------------------------------


cols_req <- c("x1", "x2", "category", "response")

m_gcm <- list()
m_gcm$name <- "gcm"
m_gcm$f_likelihood <- gcm_likelihood_no_forgetting
m_gcm$params <- params_fin$tf
m_gcm$n_feat <- 2
m_gcm$d_measure <- 1
m_gcm$lo <- lo[1:3]
m_gcm$hi <- hi[1:3]

if (is_fitting) {
  l_tbl_important_down_0 <- list()
  t_start <- Sys.time()
  l_tbl_important_down_0 <- importance_downsampling(
    l_tbl_x_ii$train[, cols_req], m_gcm, 
    # should downsampling predict responses or true category labels?
    l_tbl_x_ii$train %>% mutate(response = category),
    cat_down = 0, n_keep_max = 7#n_unique_per_category
  )
  t_end <- Sys.time()
  round(t_end - t_start, 1)
  
  saveRDS(l_tbl_important_down_0, file = "data/downsampling-ii-uniform-cat-0-third-dim.RDS")
  
  l_tbl_important_down_1 <- list()
  t_start <- Sys.time()
  l_tbl_important_down_1 <- importance_downsampling(
    l_tbl_x_ii$train[, cols_req], m_gcm, 
    # should downsampling predict responses or true category labels?
    l_tbl_x_ii$train %>% mutate(response = category),
    cat_down = 1, n_keep_max = 7#n_unique_per_category
  )
  t_end <- Sys.time()
  round(t_end - t_start, 1)
  
  saveRDS(l_tbl_important_down_1, file = "data/downsampling-ii-uniform-cat-1-third-dim.RDS")
  
  
} else {
  l_tbl_important_down_0 <- readRDS(file = "data/downsampling-ii-uniform-cat-0-third-dim.RDS")
  l_tbl_important_down_1 <- readRDS(file = "data/downsampling-ii-uniform-cat-1-third-dim.RDS")
}

extract_n_most_important_points <- function(n, l_0, l_1) {
  tbl_0 <- l_0[[n]] %>% filter(category == 0)
  tbl_1 <- l_1[[n]] %>% filter(category == 1)
  tbl_important <- rbind(tbl_0, tbl_1)
  tbl_important$trial_id <- sample(1:nrow(tbl_important), nrow(tbl_important), replace = FALSE)
  return(tbl_important)
}

l_tbl_important <- map(
  1:length(l_tbl_important_down_0), 
  l_0 = l_tbl_important_down_0, 
  l_1 = l_tbl_important_down_1,
  extract_n_most_important_points
)


ggplot(l_tbl_important[[6]], aes(x1, x2)) +
  geom_label(aes(label = rank_importance - 1, size = rank_importance)) +
  geom_abline() +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, 8)) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "") + 
  theme(strip.background = element_rect(fill = "white")) + 
  scale_color_manual(values = c("skyblue2", "tomato4"), name = "") +
  scale_size_continuous(guide = "none", range = c(2.5, 6))



# compare model predictions -----------------------------------------------


# gcm with all data points to importance models
tbl_fully_crossed <- crossing(x1 = seq(0, 7.5, by = .25), x2 = seq(0, 7.5, by = .25)) %>%
  mutate(
    category = factor(as.numeric(x1 < x2)),
    # to evaluate probability of a correct response
    response = category
  )
tbl_fully_crossed$trial_id <- sample(1:nrow(tbl_fully_crossed), nrow(tbl_fully_crossed), replace = FALSE)

preds_gcm <- category_probs(params_fin$tf, tbl_fully_crossed, l_tbl_x_ii$train, n_feat, d_measure, lo[1:3], hi[1:3])$prob_correct


wrap_category_probs <- function(x, params, tbl_transfer, tbl_train, n_feat, d_measure, lo, hi){
  cat_preds <- category_probs(params, tbl_transfer, x, n_feat, d_measure, lo, hi)
  return(cat_preds$prob_correct)
}
l_preds_gcm_downsample <- map(
  l_tbl_important, 
  wrap_category_probs, 
  params = params_fin$tf, tbl_transfer = tbl_fully_crossed, 
  n_feat = n_feat, d_measure = d_measure, lo = lo[1:3], hi = hi[1:3]
)

tbl_preds_gcm_downsample <- reduce(l_preds_gcm_downsample, cbind) %>%
  as.data.frame()
colnames(tbl_preds_gcm_downsample) <- str_c("keep = ", 1:7)


tbl_fully_crossed$pred_gcm <- preds_gcm
tbl_fully_crossed <- cbind(tbl_fully_crossed, tbl_preds_gcm_downsample)

tbl_fully_crossed %>%
  dplyr::select(x1, x2, c("pred_gcm", str_c("keep = ", 1:7))) %>%
  pivot_longer(cols = c("pred_gcm", str_c("keep = ", 1:7))) %>%
  ggplot(aes(x1, x2)) +
  geom_tile(aes(fill = value)) +
  geom_abline() +
  facet_wrap(~ name, ncol = 4) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression(x[1]), y = expression(x[2])) + 
  theme(strip.background = element_rect(fill = "white")) + 
  scale_fill_viridis_c(name = "Prop.\nCorrect")

m_cor <- cor(tbl_fully_crossed[, c("pred_gcm", str_c("keep = ", 1:7))])

tbl_cors <- tibble(
  model_1 = "gcm",
  model_2 = colnames(m_cor),
  cor = m_cor[1, ]
) %>% mutate(
  model_2 = fct_inorder(model_2, ordered = TRUE),
  model_2 = fct_relevel(model_2, "pred_gcm", after = Inf),
  model_2 = fct_recode(model_2, "GCM Base" = "pred_gcm")
)

ggplot(tbl_cors %>% filter(model_2 %in% c(str_c("keep = ", 1:7), "GCM Forget", "GCM Base")), aes(model_2, cor, group = 1)) +
  geom_line() +
  geom_point(color = "white", size = 3) +
  geom_point() + 
  theme_bw() +
  scale_x_discrete(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = "Model", y = "cor (GCM Base)") +
  theme(
    strip.background = element_rect(fill = "white"),
    text = element_text(size = 16),
    axis.text.x = element_text(angle = 90)
  ) +
  coord_cartesian(ylim = c(.25, 1))


# Category Learning Model Recovery ----------------------------------------


# recovery of strategic sampling models

c <- c(1, 2) #seq(.5, 2, length.out = 3) # 4
w <- c(.5, .75) #seq(.2, .8, length.out = 3) # 4
bias <- c(.5, .75) #seq(.2, .8, length.out = 3) # 3
n_reps <- c(2, 7)
k <- seq(2, 7, by = 1)
tbl_params <- crossing(c, w, bias, n_reps, k)

list_params <- pmap(tbl_params[, c("c", "w", "bias")], ~ list(
  not_tf = c(c = ..1, w = ..2, bias = ..3),
  tf = upper_and_lower_bounds(c(c = ..1, w = ..2, bias = ..3), lo[1:3], hi[1:3])
))
l_n_reps <- map(tbl_params$n_reps, 1)
l_ks <- map(tbl_params$k, 1)

l_info <- list(
  n_feat = n_feat, 
  d_measure = d_measure, 
  lo = lo, 
  hi = hi
)

future::plan(multisession, workers = future::availableCores() - 3)
l_results_strategic <- future_pmap(
  .l = list(l_n_reps, list_params, l_ks),
  .f = generate_and_fit, 
  tbl_train_orig = l_tbl_x_ii$train, 
  l_tbl_train_strat = l_tbl_important, 
  tbl_transfer = l_tbl_x_ii$transfer, 
  l_info = l_info,
  is_strategic = TRUE,
  .progress = TRUE,
  .options = furrr_options(seed = NULL)
)
future::plan("default")
saveRDS(l_results_strategic, file = "data/recover-strat-sampling--only-down.RDS")



# recovery of models using all data 
# i.e., gcm vanilla and gcm forgetful

c <- c(1, 2) #seq(.5, 2, length.out = 3) # 4
w <- c(.5, .75) #seq(.2, .8, length.out = 3) # 4
bias <- c(.5, .75) #seq(.2, .8, length.out = 3) # 3
delta <- c(0.00001, .9)
n_reps <- c(2, 7)
k <- seq(2, 7, by = 1)
tbl_params <- crossing(c, w, bias, delta, n_reps, k)

list_params <- pmap(tbl_params[, c("c", "w", "bias", "delta")], ~ list(
  not_tf = c(c = ..1, w = ..2, bias = ..3, delta = ..4),
  tf = upper_and_lower_bounds(c(c = ..1, w = ..2, bias = ..3, delta = ..4), lo, hi)
))
l_n_reps <- map(tbl_params$n_reps, 1)
l_ks <- map(tbl_params$k, 1)


future::plan(multisession, workers = future::availableCores() - 3)
l_results_vanilla_forgetful <- future_pmap(
  .l = list(l_n_reps, list_params, l_ks),
  .f = generate_and_fit, 
  tbl_train_orig = l_tbl_x_ii$train, 
  l_tbl_train_strat = l_tbl_important, 
  tbl_transfer = l_tbl_x_ii$transfer, 
  l_info = l_info,
  is_strategic = FALSE,
  .progress = TRUE,
  .options = furrr_options(seed = NULL)
)
future::plan("default")
saveRDS(l_results_vanilla_forgetful, file = "data/recover-vanilla-and-forgetful--only-down.RDS")





# Recognition Model Recovery ----------------------------------------------


# recognition model assumes that responses come from a mixture distribution
# memory and guessing
# important and not important representations can differ with regards to two parameters
# the probability that a response comes from memory can differ
# if responses come from memory, the precision parameter can differ











