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
library(smotefamily)
library(furrr)
library(MASS)
library(rutils)
library(e1071)


path_load <- c("utils/utils.R")
walk(path_load, source)

is_fitting <- TRUE


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

n_feat <- 2
d_measure <- 1
x1 <- seq(1, 6, by = .5)
x2 <- seq(1, 6, by = .5)
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
n_trials_total <- nrow(tbl_x_ii) * 4 # second multiplicand has to be an even number to split up train and test into equally sized data sets
tbl_samples_ii <- sample_from_grid(tbl_x_ii, n_trials_total, my_sd = .04)
idx_train <- sample(1:n_trials_total, n_trials_total/2)
tbl_train <- tbl_samples_ii[idx_train, ]
tbl_transfer <- tbl_samples_ii[!(1:n_trials_total %in% idx_train), ]
tbl_train$trial_id <- sample(1:nrow(tbl_train), nrow(tbl_train), replace = FALSE)
tbl_transfer$trial_id <- sample(1:nrow(tbl_transfer), nrow(tbl_transfer), replace = FALSE)

plot_grid(tbl_train) + geom_abline()

plot_grid(tbl_transfer) + geom_abline()

# simulate choices with simplistic difficulty function
tbl_train <- simulate_responses(tbl_train)
tbl_transfer <- simulate_responses(tbl_transfer)

if (is_fitting) {
  saveRDS(tbl_train, file = "data/tbl_train.RDS")
  saveRDS(tbl_transfer, file = "data/tbl_transfer.RDS")
} else if (!is_fitting) {
  tbl_train <- readRDS(file = "data/tbl_train.RDS")
  tbl_transfer <- readRDS(file = "data/tbl_transfer.RDS")
}

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



# fit forgetful gcm on simulated cl training choices ----------------------


params <- c(c = 1, w = .5, bias = .5, delta = .5)


params_forget <- pmap(list(params, lo, hi), upper_and_lower_bounds)
params_init_forget <- params


# takes approx 2 min to run on laptop with 300 train and transfer stimuli
# takes approx 15 seconds to run on laptop with 105 train and transfer stimuli
t_start <- Sys.time()
results_forget <- optim(
  params_init_forget,
  gcm_likelihood_forgetting,
  tbl_transfer = tbl_train,
  tbl_x = tbl_train, 
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



# save gcms, because they take some time to be fit

if (is_fitting) {
  saveRDS(params_fin, file = "data/params_fin_gcm.RDS")
  saveRDS(params_fin_forget, file = "data/params_fin_gcm_forget.RDS")
} else if (!is_fitting) {
  params_fin <- readRDS(file = "data/params_fin_gcm.RDS")
  params_fin_forget <- readRDS(file = "data/params_fin_gcm_forget.RDS")
}




# fit svm on simulated cl training choices --------------------------------


m_svm <- svm(
  response ~ x1 + x2, 
  data = tbl_train, 
  type = "C-classification", 
  kernel = "linear", 
  probability = TRUE
)

summary(m_svm)

pred_cat <- predict(m_svm, tbl_train, probability = TRUE)
predict(m_svm, tbl_train, probability = TRUE)



# fit prototype model on simulated cl training choices --------------------

m_nb <- naiveBayes(response ~ x1 + x2, data = tbl_train)








# compare predictions of gcm and svm --------------------------------------



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
    tbl_train[, cols_req], m_gcm, 
    tbl_train %>% mutate(response = category),
    cat_down = 0, n_keep_max = 10#n_unique_per_category
  )
  t_end <- Sys.time()
  round(t_end - t_start, 1)
  
  saveRDS(l_tbl_important_down_0, file = "data/downsampling-ii-uniform-cat-0.RDS")
  
  l_tbl_important_down_1 <- list()
  t_start <- Sys.time()
  l_tbl_important_down_1 <- importance_downsampling(
    tbl_train[, cols_req], m_gcm, 
    tbl_train %>% mutate(response = category),
    cat_down = 1, n_keep_max = 10#n_unique_per_category
  )
  t_end <- Sys.time()
  round(t_end - t_start, 1)
  
  saveRDS(l_tbl_important_down_1, file = "data/downsampling-ii-uniform-cat-1.RDS")
  
  
} else {
  l_tbl_important_down_0 <- readRDS(file = "data/downsampling-ii-uniform-cat-0.RDS")
  l_tbl_important_down_1 <- readRDS(file = "data/downsampling-ii-uniform-cat-1.RDS")
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
  geom_label(aes(label = rank_importance - 1)) +
  geom_abline() +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, 8))


# compare predictions -----------------------------------------------------



tbl_fully_crossed <- crossing(x1 = seq(0, 6, by = .25), x2 = seq(0, 6, by = .25)) %>%
  mutate(
    category = factor(as.numeric(x1 < x2)),
    # to evaluate probability of a correct response
    response = category
  )
tbl_fully_crossed$trial_id <- sample(1:nrow(tbl_fully_crossed), nrow(tbl_fully_crossed), replace = FALSE)



preds_svm <- attr(predict(m_svm, tbl_fully_crossed, probability = TRUE), "probabilities") %>% 
  as.data.frame() %>%
  relocate("0", .before = "1")
preds_svm$category <- tbl_fully_crossed$category
preds_svm$prob_correct <- pmap_dbl(preds_svm, ~ c(..1, ..2)[..3])
preds_gcm <- category_probs(params_fin$tf, tbl_fully_crossed, tbl_train, n_feat, d_measure, lo[1:3], hi[1:3])$prob_correct
preds_gcm_downsample <- category_probs(params_fin$tf, tbl_fully_crossed, l_tbl_important[[2]], n_feat, d_measure, lo[1:3], hi[1:3])$prob_correct
preds_gcm_forget <- category_probs(params_fin_forget$tf, tbl_fully_crossed, tbl_train, n_feat, d_measure, lo, hi)$prob_correct


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
colnames(tbl_preds_gcm_downsample) <- str_c("keep = ", 1:10)

cor(cbind(preds_gcm, tbl_preds_gcm_downsample))


preds_nb <- predict(m_nb, tbl_fully_crossed, type = "raw") %>% 
  as.data.frame() %>%
  relocate("0", .before = "1")
preds_nb$category <- tbl_fully_crossed$category
preds_nb$prob_correct <- pmap_dbl(preds_nb, ~ c(..1, ..2)[..3])

tbl_fully_crossed$pred_svm <- preds_svm$prob_correct
tbl_fully_crossed$pred_nb <- preds_nb$prob_correct
tbl_fully_crossed$pred_gcm <- preds_gcm
tbl_fully_crossed$pred_gcm_forgetful <- preds_gcm_forget
tbl_fully_crossed <- cbind(tbl_fully_crossed, tbl_preds_gcm_downsample)

tbl_fully_crossed %>%
  dplyr::select(x1, x2, str_c("keep = ", 1:10)) %>%
  pivot_longer(cols = str_c("keep = ", 1:10)) %>%
  ggplot(aes(x1, x2)) +
  geom_tile(aes(fill = value)) +
  facet_wrap(~ name, ncol = 2)


ggplot(tbl_fully_crossed %>% pivot_longer(cols = c(pred_svm, pred_gcm, pred_gcm_forgetful, pred_nb, "keep = 1")), aes(x1, x2)) +
  geom_tile(aes(fill = value)) +
  geom_abline() +
  facet_wrap(~ name) +
  scale_fill_gradient2(low = "red", high = "forestgreen", mid = "white")



m_cor <- cor(tbl_fully_crossed[, c("pred_svm", "pred_gcm", "pred_gcm_forgetful", "pred_nb", str_c("keep = ", 1:10))])

tbl_cors <- tibble(
  model_1 = "gcm",
  model_2 = colnames(m_cor),
  cor = m_cor[2, ]
) %>% mutate(
  model_2 = fct_inorder(model_2, ordered = TRUE),
  model_2 = fct_relevel(model_2, "pred_gcm_forgetful", after = Inf),
  model_2 = fct_relevel(model_2, "pred_gcm", after = Inf),
  model_2 = fct_relevel(model_2, "pred_svm", after = Inf),
  model_2 = fct_relevel(model_2, "pred_nb", after = Inf),
  model_2 = fct_recode(model_2, "GCM Forget" = "pred_gcm_forgetful"),
  model_2 = fct_recode(model_2, "GCM Base" = "pred_gcm"),
  model_2 = fct_recode(model_2, "Rule" = "pred_svm"),
  model_2 = fct_recode(model_2, "Prototype" = "pred_nb")
)

ggplot(tbl_cors %>% filter(model_2 %in% c(str_c("keep = ", 1:10), "GCM Forget", "GCM Base")), aes(model_2, cor, group = 1)) +
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
  coord_cartesian(ylim = c(.5, 1))



# downsampling svm (outdated) ---------------------------------------------


# Fitting SVM to the Training set 

m_svm <- list()
m_svm$name <- "svm"
m_svm$model <- m_svm

importance_downsampling(
  tbl_train[, cols_req], m_svm,
  tbl_train %>% mutate(response = category),
  cat_down = 0, n_keep_max = 10
)

m_svm <- svm(
  response ~ x1 + x2, 
  data = tbl_train, 
  type = "C-classification", 
  kernel = "linear", 
  probability = TRUE
)

summary(m_svm)

pred_cat <- predict(m_svm, tbl_train, probability = TRUE)
predict(m_svm, tbl_train, probability = TRUE)

tbl_fully_crossed <- crossing(x1 = seq(0, 6, by = .25), x2 = seq(0, 6, by = .25))
preds_svm <- attr(predict(m_svm, tbl_fully_crossed, probability = TRUE), "probabilities")
tbl_fully_crossed$pred_svm <- preds_svm[, 1]

ggplot(tbl_fully_crossed, aes(x1, x2)) +
  geom_point(aes(color = pred_svm)) +
  geom_label(aes(label = round(pred_svm, 2))) +
  geom_abline()

# plot difficulty to categorize item given distance from boundary
tbl_fully_crossed %>%
  mutate(
    
  )







