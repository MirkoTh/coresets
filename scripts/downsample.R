# Load packages and utils -------------------------------------------------
devtools::install_github("r-lib/conflicted")

rm(list = ls())
set.seed(43995)

library(conflicted)
library(tidyverse)
library(grid)
library(gridExtra)
library(smotefamily)
library(furrr)
library(MASS)
library(rutils)
library(docstring)

path_load <- c("utils/utils.R")
walk(path_load, source)

# training set
# transfer set = training set, because people extract important points from observed data
# downsampling set = training set
# set of parameters, possibly fitted


# idea of experiment: cl training with feedback, item recognition (half old, half new), cl transfer no feedback

# idea of script
# generate samples from two categories (category structure: inf. int.)
# half of the data is used for cl training, half for cl transfer and for item recognition
# fit plain-vanilla gcm on training responses and keep parameters
# only keep K-most important data points per category (vary K from 1 to 7)

# choice model recovery on on transfer set (downsampling training set, full training set, full training set with forgetting)
# rt model recovery on recognition responses (presented data plus additionally important points)

# model comparison:
# gcm has to be fitted on transfer set, train set can be seen as part of the model
# why? it needs to compute the pairwise similarities between all observations.
# which is not possible in the train set
# svm can be fitted on transfer set as well (assuming that performance stabilized after training)
# then, we do not really have a hold-out test set to compare the models
# but we can just use bic/aic/BF on the transfer set, that should work






# generate information integration data set -------------------------------


x1 <- seq(1, 6, by = 1)
x2 <- seq(1, 6, by = 1)
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

plot_grid(tbl_x_ii) + geom_abline()

# increase n distinct data points by sampling from
# bivariate normal using tbl_x_ii x values as means
n_trials_total <- 210
tbl_samples_ii <- sample_from_grid(tbl_x_ii, n_trials_total)
idx_train <- sample(1:n_trials_total, n_trials_total/2)
tbl_train <- tbl_samples_ii[idx_train, ]
tbl_transfer <- tbl_samples_ii[!(1:n_trials_total %in% idx_train), ]

plot_grid(tbl_train) + geom_abline()

plot_grid(tbl_transfer) + geom_abline()

# simulate choices with simplistic difficulty function
tbl_train <- simulate_responses(tbl_train)
tbl_transfer <- simulate_responses(tbl_transfer)

n_bins <- 5
tbl_train %>% 
  mutate(x1_cut = cut(x1, n_bins), x2_cut = cut(x2, n_bins)) %>%
  group_by(x1_cut, x2_cut) %>%
  summarize(prop_correct = mean(accuracy)) %>%
  ungroup() %>%
  ggplot(aes(x1_cut, x2_cut)) +
  geom_tile(aes(fill = prop_correct)) +
  geom_text(aes(label = round(prop_correct, 2))) +
  scale_fill_gradient(name = "Accuracy", low = "white", high = "forestgreen", limits = c(0, 1)) + 
  geom_abline() + ggtitle("Simulated Response Accuracy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# stimuli in Ashby, Ell, & Waldron (2003): length and orientation
# length in pixels: V1 pixels
# orientation in radians: V2 * (pi / 600)
MASS::mvrnorm(200, c(300, 400), matrix(c(8000, 7800, 7800, 8000), nrow = 2)) %>%
  as.data.frame() %>% as_tibble() %>% mutate(category = "1") %>%
  rbind(
    MASS::mvrnorm(200, c(400, 300), matrix(c(8000, 7800, 7800, 8000), nrow = 2)) %>%
  as.data.frame() %>% as_tibble() %>% mutate(category = "2") 
  ) %>%
  ggplot(aes(V1, V2, aes(group = category))) +
  geom_point(aes(shape = category)) +
  geom_abline() +
  coord_cartesian(xlim = c(0, 700), ylim = c(0, 700)) +
  theme_bw()




# fit plain-vanilla gcm on simulated cl training choices ------------------


params <- c(c = 1, w = .5, bias = .5, delta = .5)

lo <- c(0, .01, .0001, 1e-10)
hi <- c(10, .99, .9999, .999999999)
params <- pmap(list(params[1:3], lo[1:3], hi[1:3]), upper_and_lower_bounds)
params_init <- params


# takes approx 2 min to run on laptop with 300 train and transfer stimuli
# takes approx 15 seconds to run on laptop with 105 train and transfer stimuli
t_start <- Sys.time()
results_pre <- optim(
  params_init,
  gcm_likelihood_no_forgetting,
  tbl_transfer = tbl_train,
  tbl_x = tbl_train, 
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



# throw out less important points up to K ---------------------------------


cols_req <- c("x1", "x2", "category", "response")
tbl_transfer_backup <- tbl_transfer
tbl_transfer <- tbl_train %>% mutate(response = category)


m_gcm <- list()
m_gcm$name <- "gcm"
m_gcm$f_likelihood <- gcm_likelihood_no_forgetting
m_gcm$params <- params_fin$tf
m_gcm$n_feat <- 2
m_gcm$d_measure <- 1
m_gcm$lo <- lo[1:3]
m_gcm$hi <- hi[1:3]

l_tbl_important_down <- list()
t_start <- Sys.time()
l_tbl_important_down <- importance_downsampling(
  tbl_train[, cols_req], m_gcm, 
  tbl_train %>% mutate(response = category),
  cat_down = 0, n_keep_max = 10#n_unique_per_category
)
t_end <- Sys.time()
round(t_end - t_start, 1)

ggplot(tbl_drop, aes(x1, x2)) +
  geom_label(aes(label = rank_importance)) +
  geom_abline() +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 7))
ggplot(l_tbl_important_down[[10]], aes(x1, x2)) +
  geom_label(aes(label = rank_importance)) +
  geom_abline() +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 7))




# Fitting SVM to the Training set 
install.packages('e1071') 
library(e1071)

m_svm <- list()
m_svm$name <- "svm"
m_svm$model <- m_svm

importance_downsampling(
  tbl_train[, cols_req], m_svm,
  tbl_train %>% mutate(response = category),
  cat_down = 0, n_keep_max = 10
)

m_svm <- svm(
  category ~ x1 + x2, 
  data = tbl_train, 
  type = "C-classification", 
  kernel = "linear", 
  probability = TRUE
  )

summary(m_svm)

pred_cat <- predict(m_svm, tbl_train)
predict(m_svm, tbl_train, probability = TRUE)

tbl_fully_crossed <- crossing(x1 = seq(0, 6, by = .1), x2 = seq(0, 6, by = .1))
tbl_fully_crossed$pred_svm <- predict(m_svm, tbl_fully_crossed)
ggplot(tbl_fully_crossed, aes(x1, x2)) +
  geom_point(aes(color = pred_svm))



