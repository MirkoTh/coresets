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
l_x1 <- list(train = seq(1, 7, by = 1), test = seq(.5, 7.5, by = 1))
l_x2 <- list(train = seq(.5, 7.5, by = 1), test = seq(1, 7, by = 1))

train_and_test_sets <- function(x1, x2) {
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
l_tbl_x_ii <- map2(l_x1, l_x2, train_and_test_sets) %>%
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
    l_tbl_x_ii$train %>% mutate(response = category),
    cat_down = 1, n_keep_max = 10#n_unique_per_category
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




# Category Learning Model Recovery ----------------------------------------









# Recognition Model Recovery ----------------------------------------------


# recognition model assumes that responses come from a mixture distribution
# memory and guessing
# important and not important representations can differ with regards to two parameters
# the probability that a response comes from memory can differ
# if responses come from memory, the precision parameter can differ











