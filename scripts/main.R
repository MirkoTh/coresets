library(tidyverse)
library(grid)
library(gridExtra)
library(smotefamily)
library(furrr)
library(MASS)

path_load <- c("utils/utils.R")
walk(path_load, source)

# create two category structures ------------------------------------------


# squared structure as in rep change project
# information integration condition using the identity function as decision boundary

# squares
x1 <- seq(1, 10, by = 1)
x2 <- seq(1, 10, by = 1)
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


pl_sq <- plot_grid(tbl_x_sq) + geom_hline(yintercept = 5.5) + geom_vline(xintercept = 5.5)

grid.draw(arrangeGrob(pl_sq, pl_ii, nrow = 1))


# select a few training examples ------------------------------------------


# create categories with few/many examples
tbl_x_ii_inb <- tbl_x_ii %>% 
  filter(
    category %in% c(0) |
      (category %in% c(1) & x1 %in% c(1, 9) & x2 %in% c(2, 10))
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



# Upsampling --------------------------------------------------------------


# borderline smote
m_smpl_bl <- smotefamily::BLSMOTE(
  tbl_samples_ii[, c("x1", "x2")], 
  as_vector(tbl_samples_ii$category),
  K = 2
)

tbl_bl <- m_smpl_bl$syn_data
colnames(tbl_bl) <- colnames(tbl_x_sq)[1:3]
tbl_bl <- tbl_samples_ii %>% rbind(tbl_bl)
pl_bl <- plot_grid(tbl_bl) + geom_abline() + ggtitle("Borderline Sampling")


# smote
smote_k <- function(k) {
  m_smpl_smote <- smotefamily::SMOTE(
    tbl_samples_ii[, c("x1", "x2")], 
    as_vector(tbl_samples_ii$category),
    K = k
  )
  
  tbl_smote <- m_smpl_smote$syn_data
  colnames(tbl_smote) <- colnames(tbl_x_sq)[1:3]
  tbl_smote <- tbl_samples_ii %>% rbind(tbl_smote)
  plot_grid(tbl_smote) + geom_abline() + 
    ggtitle(str_c("Upsampling, K = ", k))
}


l_smote <- map(c(5, 10, 20), smote_k)

grid.draw(arrangeGrob(pl_bl, l_smote[[1]], l_smote[[2]], l_smote[[3]], nrow = 2))
  


# Base Models -------------------------------------------------------------


tbl_inb <- tbl_x_ii %>% 
  filter(
    category %in% c(0) |
      (category %in% c(1) & x1 %in% seq(1, 9, by = 4) & x2 %in% seq(2, 10, by = 4))
  )
pl_inb <- plot_grid(tbl_inb %>% mutate(category = factor(category))) + ggtitle("Stimuli") + geom_abline()


# the base models use all presented data points to predict on the transfer set

# create 10 x 10 grid of data points as a transfer set
xs <- seq(1, 10, by = 1)
tbl_transfer <- crossing(x1 = xs, x2 = xs) %>% mutate(category = fct_rev(factor(x1 > x2, labels = c(1, 0))))
l_transfer_x <- split(tbl_transfer[, c("x1", "x2")], 1:nrow(tbl_transfer))
plot_grid(tbl_transfer) + geom_abline()


# base model with bias term = .5
l_category_probs.5 <- map(l_transfer_x, gcm_base, tbl_x = tbl_inb, n_feat = 2, c = 1, w = .5, bias = .5, delta = 0, d_measure = 1)
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
l_category_probs.1 <- map(l_transfer_x, gcm_base, tbl_x = tbl_inb, n_feat = 2, c = 1, w = .5, bias = .1, delta = 0, d_measure = 1)
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
  geom_abline() + coord_cartesian(xlim = c(1, 10), ylim = c(1, 10))

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

# example using imbalanced data set sequentially presenting data



# upsample only the minority category
categories <- c(0, 1)
ns_upsample <- c(0, 100)

l_tbl_upsample_inb <- map2(categories, ns_upsample, upsample, tbl_x = tbl_inb)
tbl_inb_upsample <- l_tbl_upsample_inb %>% reduce(rbind) %>% 
  group_by(x1, x2, category) %>%
  count() %>% dplyr::select(-n) %>%
  ungroup() %>%
  mutate(
    trial_id = seq(1, nrow(.))
  )

pl_inb <- plot_grid(tbl_inb %>% mutate(category = factor(category))) + ggtitle("Stimuli") + geom_abline()
pl_inb_up <- plot_grid(tbl_inb_upsample %>% mutate(category = factor(category))) + ggtitle("Upsampled")
grid.draw(arrangeGrob(pl_inb, pl_inb_up, nrow = 1))


gcm_base(tibble(x1 = 2.5, x2 = 4), tbl_inb, 2, c = 1, w = .5, bias = 1/3, delta = 0, d_measure = 1)
gcm_base(tibble(x1 = 2.5, x2 = 4), tbl_inb_upsample, 2, c = 1, w = .5, bias = 1/3, delta = 0, d_measure = 1)



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


# for one category
categories <- c(0, 1)
ns_upsample <- c(0, 100)

l_tbl_upsample_inb <- map2(categories, ns_upsample, upsample, tbl_x = tbl_inb)
# add small amount of noise to circumvent ties in rank odering
tbl_inb_upsample <- l_tbl_upsample_inb %>% reduce(rbind) %>% 
  group_by(x1, x2, category) %>%
  count() %>% dplyr::select(-n) %>%
  ungroup() %>%
  mutate(
    trial_id = seq(1, nrow(.)),
    x1 = x1 + rnorm(nrow(.), 0, .01),
    x2 = x2 + rnorm(nrow(.), 0, .01)
  )
plot_grid(tbl_inb_upsample %>% mutate(category = factor(category))) + geom_abline()
params <- c(c = 1, w = .5, bias = .5)


# calculate importance given set of data points
tbl_inb_plus <- tbl_inb %>% dplyr::select(x1, x2, category)
tbl_importance <- tbl_inb_upsample
l_new_samples <- split(tbl_inb_upsample %>% dplyr::select(-trial_id), tbl_inb_upsample$trial_id)

tbl_important_samples_2 <- importance_upsampling(l_new_samples, tbl_importance, tbl_inb_plus, params, tbl_transfer, n_feat = 2, d_measure = 1, n_max = 2)
tbl_important_samples_2$is_new <- 1
tbl_important_samples_2$is_new[(nrow(tbl_inb) + 1) : nrow(tbl_important_samples_2)] <- 5

tbl_important_samples_10 <- importance_upsampling(l_new_samples, tbl_importance, tbl_inb_plus, params, tbl_transfer, n_feat = 2, d_measure = 1, n_max = 10)
tbl_important_samples_10$is_new <- 1
tbl_important_samples_10$is_new[(nrow(tbl_inb) + 1) : nrow(tbl_important_samples_10)] <- 5


pl_points_obs <- plot_grid(tbl_inb %>% mutate(category = factor(category))) + 
  geom_abline() + ggtitle("Observed")
pl_points_transfer <- plot_grid(tbl_transfer) + geom_abline() + ggtitle("Transfer")
pl_points_importance_2 <- plot_grid(tbl_important_samples_2) + geom_abline() +
  geom_point(aes(size = is_new), shape = 1) +
  scale_size_continuous(guide = "none") +
  ggtitle("Add Two Samples")
pl_points_importance_10 <- plot_grid(tbl_important_samples_10) + geom_abline() +
  geom_point(aes(size = is_new), shape = 1) +
  scale_size_continuous(guide = "none") +
  ggtitle("Add Ten Samples")

grid.draw(arrangeGrob(pl_points_obs, pl_points_transfer, pl_points_importance_2, pl_points_importance_10, nrow = 2))




plot_grid(tbl_drop)

# todos
# up- and downsampling only for relevant category (i.e., lo freq and hi freq, respectively)
# assumption: up- and downsampling are achieved using originally presented stimuli




