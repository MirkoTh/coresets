rm(list = ls())
set.seed(43995)

library(tidyverse)
library(grid)
library(gridExtra)
library(smotefamily)
library(furrr)
library(MASS)
library(docstring)

path_load <- c("utils/utils.R")
walk(path_load, source)

# create two category structures ------------------------------------------


# squared structure as in rep change project
# information integration condition using the identity function as decision boundary

# squares
x1 <- seq(1, 6, by = 1)
x2 <- seq(1, 6, by = 1)
tbl_x_ii <- crossing(x1, x2)
tbl_x_ii$category <- factor(as.numeric(tbl_x_ii$x1 < tbl_x_ii$x2))
tbl_x_ii <- tbl_x_ii[!(tbl_x_ii$x1 == tbl_x_ii$x2), ]
tbl_x_ii$cat_structure <- "information-integration"
tbl_x_ii$trial_id <- seq(1, nrow(tbl_x_ii), by = 1)

pl_ii <- plot_grid(tbl_x_ii) + geom_abline()



# square
tbl_x_sq <- crossing(x1, x2)
tbl_x_sq$category <- factor((tbl_x_sq$x1 > mean(x1)) + (tbl_x_sq$x2 > mean(x2)) * 2 )
tbl_x_sq$cat_structure <- "square"
tbl_x_sq$trial_id <- seq(1, nrow(tbl_x_sq), by = 1)


pl_sq <- plot_grid(tbl_x_sq) + geom_hline(yintercept = mean(tbl_x_sq$x1)) + geom_vline(xintercept = mean(tbl_x_sq$x2))

grid.draw(arrangeGrob(pl_sq, pl_ii, nrow = 1))


# select a few training examples ------------------------------------------


# create categories with few/many examples
tbl_x_ii_inb <- tbl_x_ii %>% 
  filter(
    category %in% c(0) |
      (category %in% c(1) & x1 %in% c(1, 5) & x2 %in% c(2, 6))
    #| (category == 1 & x1 == 8 & x2 == 3) 
    # | (category == 2 & x1 == 3 & x2 == 8) 
  ) %>% mutate(
    n = 10
  )
tbl_x_ii_inb <- tbl_x_ii_inb[
  sample(1:nrow(tbl_x_ii_inb), nrow(tbl_x_ii_inb), replace = FALSE), 
]
tbl_x_ii_inb$trial_id <- seq(1, nrow(tbl_x_ii_inb), by = 1)
# tbl_x_sq_inb$n[tbl_x_sq_inb$category %in% c(1, 2)] <- 3

plot_grid(tbl_x_ii_inb) + geom_abline()



l_samples_ii <- pmap(tbl_x_ii_inb[, c("x1", "x2", "category", "n")], sample_2d_stimuli, sd = .2)
tbl_samples_ii <- reduce(l_samples_ii, rbind)

plot_grid(tbl_samples_ii) + geom_abline()


# Base Models -------------------------------------------------------------


tbl_imb <- tbl_x_ii %>% 
  filter(
    category %in% c(0) |
      (category %in% c(1) & x1 %in% seq(1, 5, by = 4) & x2 %in% seq(2, 6, by = 4))
  ) %>% mutate(trial_id = sample(1:nrow(.), nrow(.), replace = FALSE)) %>%
  arrange(trial_id)
tbl_imb_weighted <- tbl_imb %>%
  rbind(
    repeat_tibble(
      tbl_imb %>% filter(category == 1), 4
    )
  ) %>% mutate(trial_id = sample(1:nrow(.), nrow(.), replace = FALSE)) %>%
  arrange(trial_id)


pl_inb <- plot_grid(tbl_imb_weighted %>% mutate(category = factor(category))) + ggtitle("Stimuli") + geom_abline()


# the base models use all presented data points to predict on the transfer set

# create 10 x 10 grid of data points as a transfer set

tbl_transfer <- tbl_x_ii %>%
  mutate(
    trial_id = sample(1:nrow(.), nrow(.), replace = FALSE)
  ) %>% arrange(trial_id)
l_transfer_x <- split(tbl_transfer[, c("x1", "x2")], 1:nrow(tbl_transfer))
plot_grid(tbl_transfer) + geom_abline()

# simulate response

# add a bit of noise to circumvent ties in rank ordering
# when calculating importance

tbl_imb[, c("x1", "x2", "category")] <- add_jitter(tbl_imb)
tbl_imb_weighted[, c("x1", "x2", "category")] <- add_jitter(tbl_imb_weighted)
tbl_transfer[, c("x1", "x2", "category")] <- add_jitter(tbl_transfer)


tbl_imb <- simulate_responses(tbl_imb)
tbl_imb_weighted <- simulate_responses(tbl_imb_weighted)
tbl_transfer <- simulate_responses(tbl_transfer)

# tbl_imb_old <- tbl_imb
# tbl_imb_weighted_old <- tbl_imb_weighted
# tbl_transfer_old <- tbl_transfer
# 
# tbl_transfer$response == tbl_transfer_old$response
# tbl_imb$response == tbl_imb_old$response
# tbl_imb_weighted$response == tbl_imb_weighted_old$response

# base model with bias term = .5
l_category_probs.5 <- map(l_transfer_x, gcm_base, tbl_x = tbl_imb, n_feat = 2, c = 1, w = .5, bias = .5, delta = 0, d_measure = 1)
tbl_category_probs.5 <- cbind(
  tbl_transfer, 
  as_tibble(as.data.frame(reduce(l_category_probs.5, rbind)))
) %>% as_tibble()
tbl_category_probs.5$prob_error <- pmap_dbl(
  tbl_category_probs.5[, c("0", "1", "category")],
  ~ 1 - c(..1, ..2)[as.numeric(as.character(..3)) + 1]
)
pl_gcm_baseline.5 <- plot_grid(tbl_category_probs.5) +
  geom_tile(aes(fill = prob_error), alpha = .5) +
  geom_text(aes(label = round(prob_error, 2))) +
  scale_fill_gradient(name = "Prob. Error", low = "white", high = "tomato", limits = c(0, 1)) + 
  geom_abline() + ggtitle("Baseline Classification, Bias = .5")

# base model with bias term = .1
l_category_probs.1 <- map(l_transfer_x, gcm_base, tbl_x = tbl_imb, n_feat = 2, c = 1, w = .5, bias = .1, delta = 0, d_measure = 1)
tbl_category_probs.1 <- cbind(
  tbl_transfer, 
  as_tibble(as.data.frame(reduce(l_category_probs.1, rbind)))
) %>% as_tibble()
tbl_category_probs.1$prob_error <- pmap_dbl(
  tbl_category_probs.1[, c("0", "1", "category")],
  ~ 1 - c(..1, ..2)[as.numeric(as.character(..3)) + 1]
)
pl_gcm_baseline.1 <- plot_grid(tbl_category_probs.1) +
  geom_tile(aes(fill = prob_error), alpha = .5) +
  geom_text(aes(label = round(prob_error, 2))) +
  scale_fill_gradient(name = "Prob. Error", low = "white", high = "tomato", limits = c(0, 1)) + 
  geom_abline() + ggtitle("Baseline Classification, Bias = .1")

# plot together
grid.draw(arrangeGrob(pl_inb, pl_gcm_baseline.5, pl_gcm_baseline.1, nrow = 1))



# Upsampling Example ------------------------------------------------------

x_new <- tibble(x1 = 4.5, x2 = 4.7)
tbl_x <- tbl_x_ii
n_feat <- 2
c <- .5
w <- 1/n_feat # equal
bias <- .5
delta <- .99
d_measure <- 1

gcm_base(x_new, tbl_x, 2, c, w, bias, delta, d_measure)
gcm_base(x_new, tbl_x, 2, c, w, bias, 0, d_measure)


plot_grid(tbl_x) +
  geom_abline() + coord_cartesian(xlim = c(1, max(tbl_x$x1)), ylim = c(1, max(tbl_x$x2)))

# for one category
categories <- c(0, 1)
ns_upsample <- c(0, 1000)

l_tbl_upsample <- map2(categories, ns_upsample, upsample, tbl_x = tbl_x)
tbl_upsample <- l_tbl_upsample %>% reduce(rbind) %>% mutate(trial_id = seq(1, nrow(.)))
plot_grid(tbl_upsample %>% mutate(category = factor(category)))
gcm_base(x_new, tbl_upsample, 2, c, w, 0, d_measure)




# Upsampling & Importance Sampling ----------------------------------------


# upsample minority category
# add/remove items from the set of presented items using importance sampling
# this makes sense only in the non-decay model, because importance or weight is given by the decay parameter
# number of finally used items is controlled by capacity parameter


# todo
# create grid of 10 x 10 values crossing x1 and x2
# 1. run gcm with upsampled data sets (few added points - many added points)
# 2. run gcm with bias parameter varying from unbiased to very biased and plot probabilities category 0 over 100 grid points

# main remaining issue: how is upsampling exactly achieved in the model? given I want to add n data points to one category
# possible solution: do upsampling as currently implemented in upsample function, 
# and then rank sampled data points according to their importance for the categorization problem
# importance sampling has to be done sequentially, otherwise points from the same neighbourhood will be sampled
# only keep the n most important data points
# adjust bias parameter accordingly, such that ratio of ndata(cat_0)/ndata(cat_1) is put into bias parameter
# e.g., cat0 50, cat1 100 --> bias = c(.66, .33)

# upsampling
# need parameter capacity: nr data points upsampled
# write function calculating likelihood of 10x10 data points
# when each of the upsampled data points is added to the original data set
# pick the best point and repeat process
# pick the best point and ... --> repeat nr data points upsampled

# general idea:
# what information use people to draw inferences when they have comparatively little/much information about a category?
# a) gcm: they just use the presented information as is
# b) prototype: they use the average representation of the category
# maybe just focus on hypotheses in c), because unclear from previous research whether a) or b)?
# c) they do not use all the information presented to them
# 1. because they forget according to a decay function
# 2. because they strategically remember/sample some information: up-sampling & down-sampling given a capacity limit


# varying number of upsampled points and evaluating likelihood of given responses
# possible solution: first, upsample sufficient number of data points (e.g., 50)
# then, do a constrained optimization with n ranging between 0 and that number (i.e., 50)



# first, upsample sufficient number of data points for the minority category
# for one category
categories <- c(0, 1)
ns_upsample <- c(0, 10)

l_tbl_upsample_inb <- map2(categories, ns_upsample, upsample, tbl_x = tbl_imb)
# add small amount of noise to circumvent ties in rank ordering
tbl_imb_upsample <- l_tbl_upsample_inb %>% reduce(rbind) %>% 
  mutate(category = as.factor(category)) %>%
  left_join(tbl_imb[, c("x1", "x2", "category", "accuracy")], by = c("x1", "x2", "category")) %>%
  mutate(
    is_new = is.na(accuracy)
  ) %>% dplyr::select(-accuracy) %>%
  group_by(x1, x2, category, is_new) %>%
  count() %>% dplyr::select(-n) %>%
  ungroup() %>%
  mutate(
    trial_id = seq(1, nrow(.)),
    x1 = x1 + rnorm(nrow(.), 0, .01),
    x2 = x2 + rnorm(nrow(.), 0, .01),
    response = category
  )
plot_grid(tbl_imb_upsample %>% mutate(category = factor(category))) + geom_abline()



pl_train <- plot_grid(tbl_imb)
pl_transfer <- plot_grid(tbl_transfer)
grid.draw(arrangeGrob(pl_train, pl_transfer, nrow = 1))





# params have to be chosen
# first fit the model on the presented data and then fix the parameters
# then up- or downsample


# to be tested
# does it make a difference on the specific values of the fit model
# which data points are considered to be important?
# e.g., 1, .5, .5 vs. a model fit on a given set of data (as below)

#gcm_likelihood_no_forgetting(params_init, tbl_transfer, tbl_imb, 2, 1, lo, hi)

params <- c(c = 1, w = .5, bias = .5)

lo <- c(0, .01, .0001)
hi <- c(10, .99, .9999)
params <- pmap(list(params, lo, hi), upper_and_lower_bounds)
params_init <- params



t_start <- Sys.time()
results_pre <- optim(
  params_init,
  gcm_likelihood_no_forgetting,
  tbl_transfer = tbl_imb_weighted,
  tbl_x = tbl_imb_weighted, 
  n_feat = 2,
  d_measure = 1,
  lo = lo,
  hi = hi
)
t_end <- Sys.time()
round(t_end - t_start, 1)

params_fin <- list()
params_fin[["not_tf"]] <- pmap_dbl(list(results_pre$par, lo, hi), upper_and_lower_bounds_revert)
params_fin[["tf"]] <- results_pre$par
params_ignorant <- list()
params_ignorant[["tf"]] <- pmap_dbl(list(c(c = 1, w = .5, bias = .5), lo, hi), upper_and_lower_bounds)
params_ignorant[["not_tf"]] <- c(c = 1, w = .5, bias = .5)
params_fin[["not_tf"]]
params_ignorant[["not_tf"]]



tbl_imb_weighted %>%
  mutate(x1 = round(x1, 0), x2 = round(x2, 0)) %>%
  group_by(x1, x2) %>%
  summarize(accuracy = mean(accuracy)) %>%
  ggplot(aes(x1, x2)) +
  geom_tile(aes(fill = accuracy))
# training data
# 15 vs. 3 examples, but number of presentations per category is fixed


# two options
# 1. present a balanced data set for training (with fewer examples in one category)
# 2. present an imbalanced data set for training (with fewer presentations and examples in one category)

# for option 1., the upsampled points can be matched in the "number of presentations" as compared to the minority examples
# or the upsampled points are just added as one presentation, which reduces their influence compared to 1.
# for option 2., this is not a problem. but only presenting a few stimuli in one category may be odd to start with



# First Up, Then Down -----------------------------------------------------


n_unique_per_category <- 9
tbl_n_change <- tbl_imb_weighted %>%
  mutate(x1 = round(x1, 0), x2 = round(x2, 0)) %>%
  count(x1, x2, category) %>% count(category) %>%
  mutate(n_change = n_unique_per_category - n)

cols_req <- c("x1", "x2", "category", "response")

# upsampling only on minority category
tbl_importance <- tbl_imb_upsample %>% 
  filter(category == 1 & is_new) %>% 
  dplyr::select(-is_new)
tbl_important_up <- importance_upsampling(
  tbl_importance, tbl_imb_weighted[, cols_req], params_ignorant$tf, 
  tbl_transfer %>% mutate(response = category),
  n_feat = 2, d_measure = 1, lo = lo, hi = hi, n_max = 6
)


# downsampling only on majority category
tbl_important_down <- importance_downsampling(
  tbl_important_up[, cols_req], params_fin$tf, 
  tbl_transfer %>% mutate(response = category), n_feat = 2, d_measure = 1, 
  lo = lo, hi = hi, cat_down = 0, n_max = n_unique_per_category
)

plot_grid(tbl_important_down)

tbl_all <- mark_changes(tbl_important_down, tbl_imb_weighted)


pl_points_obs <- plot_grid(tbl_imb %>% mutate(category = factor(category))) + 
  geom_abline() + ggtitle("Observed")
pl_points_transfer <- plot_grid(tbl_transfer) + geom_abline() + ggtitle("Transfer")
pl_points_importance <- plot_new_and_dropped(tbl_all)

grid.draw(arrangeGrob(pl_points_obs, pl_points_transfer, pl_points_importance, nrow = 1))
