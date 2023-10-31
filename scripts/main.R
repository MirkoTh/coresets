rm(list = ls())
set.seed(43995)

library(tidyverse)
library(grid)
library(gridExtra)
library(smotefamily)
library(furrr)
library(MASS)
library(docstring)
library(rutils)

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
tbl_x_ii_imb <- tbl_x_ii %>% 
  filter(
    category %in% c(0) |
      (category %in% c(1) & x1 %in% c(1, 5) & x2 %in% c(2, 6))
    #| (category == 1 & x1 == 8 & x2 == 3) 
    # | (category == 2 & x1 == 3 & x2 == 8) 
  )
tbl_x_ii_imb <- tbl_x_ii_imb[
  sample(1:nrow(tbl_x_ii_imb), nrow(tbl_x_ii_imb), replace = FALSE), 
]
tbl_x_ii_imb$trial_id <- seq(1, nrow(tbl_x_ii_imb), by = 1)

n_trials_total <- 300
n_stim <- tbl_x_ii_imb %>% count(category)
n_stimulus_reps <- tibble(
  category = n_stim$category,
  n = (n_trials_total / 2) / n_stim$n
)

tbl_x_ii_imb <- tbl_x_ii_imb %>% left_join(n_stimulus_reps, by = "category")



#tbl_x_sq_imb$n[tbl_x_sq_imb$category %in% c(1, 2)] <- 3

plot_grid(tbl_x_ii_imb) + geom_abline()


tbl_x_ii_imb %>% count(category)

overlap <- 1

# do not allow samples to come from wrong category
while(overlap > 0) {
  l_samples_ii <- pmap(tbl_x_ii_imb[, c("x1", "x2", "category", "n")], sample_2d_stimuli, sd = .075)
  tbl_samples_ii <- reduce(l_samples_ii, rbind)
  tbl_overlap <- tbl_samples_ii %>% group_by(category) %>%
    count(cat_wrong_1 = x1 >= x2, cat_wrong_0 = x2 >= x1) %>%
    pivot_longer(cols = c(cat_wrong_0, cat_wrong_1)) %>%
    mutate(cat_wrong = str_extract(name, "[0-9]$")) %>%
    filter(value & as.numeric(as.character(category)) == cat_wrong)
  if (nrow(tbl_overlap) == 0) {overlap <- 0}
}
tbl_samples_ii <- tbl_samples_ii %>%
  mutate(trial_id = sample(1:nrow(.), nrow(.)))


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

# use samples instead of data points on the grid
tbl_imb <- tbl_samples_ii
tbl_imb_weighted <- tbl_samples_ii


pl_imb <- plot_grid(tbl_imb_weighted %>% mutate(category = factor(category))) + ggtitle("Stimuli") + geom_abline()


# the base models use all presented data points to predict on the transfer set

# create representative grid of data points as transfer set

tbl_transfer <- tbl_x_ii %>%
  mutate(
    trial_id = sample(1:nrow(.), nrow(.), replace = FALSE)
  ) %>% arrange(trial_id)
l_transfer_x <- split(tbl_transfer[, c("x1", "x2")], 1:nrow(tbl_transfer))
plot_grid(tbl_transfer) + geom_abline()

# simulate response

# add a bit of noise to circumvent ties in rank ordering
# when calculating importance
# 
# tbl_imb[, c("x1", "x2", "category")] <- add_jitter(tbl_imb)
# tbl_imb_weighted[, c("x1", "x2", "category")] <- add_jitter(tbl_imb_weighted)
tbl_transfer[, c("x1", "x2", "category")] <- add_jitter(tbl_transfer)


tbl_imb <- simulate_responses(tbl_imb)
tbl_imb_weighted <- simulate_responses(tbl_imb_weighted)
tbl_transfer <- simulate_responses(tbl_transfer)

saveRDS(tbl_transfer, file = "data/transfer-data.RDS")
saveRDS(tbl_imb, file = "data/hotspot-data.RDS")

# The Effect of Category Bias ---------------------------------------------


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
grid.draw(arrangeGrob(pl_imb, pl_gcm_baseline.5, pl_gcm_baseline.1, nrow = 1))



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
  geom_abline() + coord_cartesian(xlim = c(1, max(tbl_x$x1)), ylim = c(1, max(tbl_x$x2))) +
  ggstar::geom_star(data = x_new %>% mutate(category = 3), aes(x1, x2))

# for one category
categories <- c(0, 1)
ns_upsample <- c(0, 1000)

l_tbl_upsample <- map2(categories, ns_upsample, upsample, tbl_x = tbl_x) # tbl_samples_ii
tbl_upsample <- l_tbl_upsample %>% reduce(rbind) %>% mutate(trial_id = seq(1, nrow(.)))
plot_grid(tbl_upsample %>% mutate(category = factor(category)))
gcm_base(x_new, tbl_upsample, 2, c, w, 0, d_measure)




# Upsampling & Importance Sampling ----------------------------------------


# upsample minority category
# add/remove items from the set of presented items using importance sampling
# this makes sense only in the non-decay model, because importance or weight is given by the decay parameter
# number of finally used items is controlled by capacity parameter


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




clusters_minority <- kmeans(tbl_imb[tbl_imb$category == 1, c("x1", "x2")], 3)
tbl_clusters <- tibble(as.data.frame(clusters_minority$centers)) %>%
  mutate(category = 1) %>%
  rbind(tbl_imb %>% filter(category == 0) %>% dplyr::select(x1, x2, category))
tbl_clusters <- simulate_responses(tbl_clusters)
plot_grid(tbl_clusters) + geom_abline()

# first, upsample sufficient number of data points for the minority category
# for one category
categories <- c(0, 1)
ns_upsample <- c(0, 9)

l_tbl_upsample_imb <- map2(categories, ns_upsample, upsample, tbl_x = tbl_clusters)
# add small amount of noise to circumvent ties in rank ordering

# left_join on tbl_imb before using tbl_clusters
tbl_imb_upsample <- l_tbl_upsample_imb %>% reduce(rbind) %>% 
  mutate(category = as.factor(category)) %>%
  left_join(tbl_clusters[, c("x1", "x2", "category", "accuracy")], by = c("x1", "x2", "category")) %>%
  mutate(
    is_new = is.na(accuracy)
  ) %>% dplyr::select(-accuracy) %>%
  group_by(x1, x2, category, is_new) %>%
  count() %>% dplyr::select(-n) %>%
  ungroup() %>%
  mutate(
    trial_id = seq(1, nrow(.)),
    # x1 = x1 + rnorm(nrow(.), 0, .01),
    # x2 = x2 + rnorm(nrow(.), 0, .01),
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

params <- c(c = 1, w = .5, bias = .5, delta = .5)

lo <- c(0, .01, .0001, 1e-10)
hi <- c(10, .99, .9999, .999999999)
params <- pmap(list(params[1:3], lo[1:3], hi[1:3]), upper_and_lower_bounds)
params_init <- params


# takes approx 1 min to run on lab computer
t_start <- Sys.time()
results_pre <- optim(
  params_init,
  gcm_likelihood_no_forgetting,
  tbl_transfer = tbl_imb_weighted,
  tbl_x = tbl_imb_weighted, 
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


n_unique_per_category <- 10
tbl_n_change <- tbl_imb_weighted %>%
  mutate(x1 = round(x1, 0), x2 = round(x2, 0)) %>%
  count(x1, x2, category) %>% count(category) %>%
  mutate(n_change = n_unique_per_category - n)

cols_req <- c("x1", "x2", "category", "response")

# upsampling only on minority category
# takes approx. 10 seconds on lab computer

tbl_importance <- tbl_imb_upsample %>% 
  filter(category == 1 & is_new) %>% 
  dplyr::select(-is_new)
n_add <- n_unique_per_category - nrow(tbl_imb_upsample %>% filter(category == 1 & !is_new))

t_start <- Sys.time()
l_tbl_important_up <- importance_upsampling(
  tbl_importance, tbl_imb_weighted[, cols_req], params_fin$tf, 
  tbl_transfer %>% mutate(response = category),
  n_feat = 2, d_measure = 1, lo = lo[1:3], hi = hi[1:3], n_add = n_add
)
t_end <- Sys.time()
round(t_end - t_start, 1)

tbl_important_up <- l_tbl_important_up[[n_add]]
plot_grid(tbl_important_up)



# downsampling only on majority category
# takes about 6.5 mins for 150 stimuli in majority category on laptop
# takes about 8.4 mins for 150 stimuli in majority category on lab computer
l_tbl_important_down <- list()
t_start <- Sys.time()
l_tbl_important_down <- importance_downsampling(
  tbl_imb_weighted[, cols_req], params_fin$tf, 
  tbl_transfer %>% mutate(response = category), n_feat = 2, d_measure = 1, 
  lo = lo[1:3], hi = hi[1:3], cat_down = 0, n_keep_max = n_unique_per_category
)
t_end <- Sys.time()
round(t_end - t_start, 1)

tbl_important_down <- l_tbl_important_down[[n_unique_per_category]]

l_tbl_up_and_down <- list()
l_tbl_changes_up <- list()
l_tbl_changes_down <- list()


for (i in 1:length(l_tbl_important_up)) {
  tbl_up_and_down <- rbind(l_tbl_important_up[[i]], l_tbl_important_down[[i + 3]] %>% dplyr::select(-c(importance, rank_importance)))
  l_tbl_changes_up[[i]] <- mark_changes(l_tbl_important_up[[i]], tbl_imb_weighted %>% mutate(cat_structure = "Information Integration"))
  l_tbl_changes_down[[i]] <- mark_changes(l_tbl_important_down[[i + 3]], tbl_imb_weighted %>% mutate(cat_structure = "Information Integration"))
  l_tbl_up_and_down[[i]] <- l_tbl_changes_up[[i]] %>% filter(is_new) %>%
    rbind(l_tbl_changes_down[[i]] %>% filter(!is_dropped & category == 0)) %>%
    rbind(tbl_clusters %>% filter(category == 1) %>% dplyr::select(x1, x2, category) %>% mutate(is_new = FALSE, is_dropped = FALSE))
  
}

saveRDS(l_tbl_up_and_down, "data/l_tbl_up_and_down.RDS")


plot_grid(tbl_up_and_down)

tbl_important_down <- tbl_important_down %>%
  mutate(trial_id = sample(1:nrow(.), nrow(.), replace = FALSE)) %>%
  arrange(trial_id)

plot_grid(tbl_important_down)

pl_points_obs <- plot_grid(tbl_imb %>% mutate(category = factor(category))) + 
  geom_abline() + ggtitle("Observed")
pl_points_transfer <- plot_grid(tbl_transfer) + geom_abline() + ggtitle("Transfer")
pl_points_importance <- plot_new_and_dropped(l_tbl_up_and_down[[7]] %>% rbind(l_tbl_changes_down[[7]] %>% filter(category == 0)))

grid.draw(arrangeGrob(pl_points_obs, pl_points_transfer, pl_points_importance, nrow = 1))




# Model Recovery ----------------------------------------------------------

# todos
# compare strategic sampling model with decay gcm model, check that trial numbers of data points in training data set are ok
# recover models with more data points (aka hot spots vs. grid with possibly 100 examples per category)




# Predict Data with Strategic Sampling Model


# Test Model Recovery
# Model Candidates:
# - strategic sampling (true)
# - using all presented data points
# - decay gcm (not yet implemented)
# diffs in lls
# 10 reps: 6
# 15 reps: 
# 20 reps: 20

# data required to put the model-recovery into a separate file:
# tbl_imb_weighted
# l_tbl_up_and_down
l_tbl_up_and_down <- readRDS("data/l_tbl_up_and_down.RDS")
tbl_transfer



plot_grid(tbl_imb_weighted)
plot_grid(tbl_up_and_down)

c <- seq(.5, 2, length.out = 4)
w <- seq(.2, .8, length.out = 4)
bias <- seq(.2, .8, length.out = 3)
n_reps <- c(5, 10, 20)
k <- seq(1, 7, by = 1)
tbl_params <- crossing(c, w, bias, n_reps, k)

list_params <- pmap(tbl_params[, c("c", "w", "bias")], ~ list(
  not_tf = c(c = ..1, w = ..2, bias = ..3),
  tf = upper_and_lower_bounds(c(c = ..1, w = ..2, bias = ..3), lo[1:3], hi[1:3])
))
l_n_reps <- map(tbl_params$n_reps, 1)
l_ks <- map(tbl_params$k, 1)

future::plan(multisession, workers = future::availableCores() - 3)
l_results <- future_pmap(
  .l = list(l_n_reps, list_params, l_ks),
  .f = generate_and_fit, 
  tbl_train_orig = tbl_imb_weighted, 
  l_tbl_train_strat = l_tbl_up_and_down, 
  tbl_transfer = tbl_transfer, 
  n_feat = n_feat, 
  d_measure = d_measure, 
  lo = lo, 
  hi = hi,
  .progress = TRUE,
  .options = furrr_options(seed = NULL)
)
saveRDS(l_results, file = "data/recovery-hotspots.RDS")
future::plan("default")


#l_results <- readRDS(file = "data/recovery-hotspots.RDS")


# comparisons
# 1. strategic in vs. strategic out, original out, and decay out bic -> can we recover the generating nr of examples
# 2. strategic gen vs. strategic out params -> can we recover the parameters


# unwrap lls from different number of strategically sampled data points
l_results_ll <- map(l_results, "n2lls")
tbl_ll_strat <- reduce(map(1:7, ~ map_dbl(map(l_results_ll, 1), .x)), cbind) %>%
  as.data.frame()
colnames(tbl_ll_strat) <- str_c("k", k)
ll_original <- map_dbl(l_results_ll, "original")
ll_decay <- map_dbl(l_results_ll, "decay")

tbl_lls <- cbind(
  tbl_params, 
  tbl_ll_strat, ll_original, ll_decay) %>%
  mutate(
    bic_k1 = 3*log(nrow(tbl_transfer)) + k1,
    bic_k2 = 3*log(nrow(tbl_transfer)) + k2,
    bic_k3 = 3*log(nrow(tbl_transfer)) + k3,
    bic_k4 = 3*log(nrow(tbl_transfer)) + k4,
    bic_k5 = 3*log(nrow(tbl_transfer)) + k5,
    bic_k6 = 3*log(nrow(tbl_transfer)) + k6,
    bic_k7 = 3*log(nrow(tbl_transfer)) + k7,
    bic_orig = 3*log(nrow(tbl_transfer)) + ll_original,
    bic_decay = 4*log(nrow(tbl_transfer)) + ll_decay
  )


tbl_bics_long <- tbl_lls %>%
  as_tibble() %>%
  rename(gen_k = k) %>%
  dplyr::select(-c(starts_with("k"), "ll_original", "ll_decay")) %>%
  pivot_longer(starts_with("bic"))
tbl_bics_long$gen_k <- factor(tbl_bics_long$gen_k, labels = str_c("k = ", 1:7), ordered = TRUE)
tbl_bics_long$name <- fct_inorder(factor(tbl_bics_long$name, ordered = TRUE))
levels(tbl_bics_long$name) <- c(str_c("k = ", 1:7), "Presented", "Decay")
tbl_bics_long <- tbl_bics_long %>% 
  group_by(c, w, bias, n_reps, gen_k) %>%
  mutate(bic_prop = value / max(value)) %>%
  ungroup()
tbl_bics_long_agg <- tbl_bics_long %>% 
  group_by(n_reps, gen_k, name) %>%
  summarize(bic_prop_mn = mean(bic_prop)) %>%
  ungroup()
tbl_winner <- tbl_bics_long_agg %>%
  group_by(gen_k, n_reps) %>%
  mutate(bic_rwn = row_number(bic_prop_mn)) %>%
  filter(bic_rwn == 1) %>%
  ungroup()

ggplot(tbl_bics_long_agg, aes(gen_k, name)) + 
  geom_tile(aes(fill = bic_prop_mn), color = "black") +
  geom_label(aes(label = round(100 * bic_prop_mn, 1))) +
  geom_tile(data = tbl_winner, aes(gen_k, name), fill = NA, color = "white", size = 2) +
  theme_bw() +
  facet_wrap(~ n_reps) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Generating Model", y = "Recovered Model", title = "100 * BIC / BIC (Worst) by Gen. Model") + 
  theme(strip.background = element_rect(fill = "white")) + 
  scale_fill_gradient(low = "skyblue2", high = "tomato4", name = "", guide = "none")



# 
# tbl_prop_success <- grouped_agg(tbl_lls %>% mutate(success = bic_delta < 0), n_reps, success)
# 
# ggplot(tbl_lls, aes(bic_delta)) + 
#   geom_histogram(color = "white", fill = "skyblue2") +
#   geom_vline(xintercept = 0, color = "black", linetype = "dotdash") +
#   geom_label(
#     data = tbl_prop_success, 
#     aes(x = -mean(tbl_lls$bic_delta), y = nrow(tbl_lls)/20, 
#         label = str_c("Prop. Recov. = ", round(mean_success, 2))
#     )
#   ) + facet_wrap(~ n_reps) +
#   theme_bw() +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x = "-2*LL Delta (Strat. Sampl. - Original)", y = "Nr. Simulations") + 
#   theme(strip.background = element_rect(fill = "white")) + 
#   scale_color_manual(values = c("skyblue2", "tomato4"), name = "")


l_results_params <- map(l_results, "params_strat")
l_results_params_gen <- map2(1:nrow(tbl_params), tbl_params$k, ~ l_results_params[[.x]][[.y]][["not_tf"]])
tbl_params_fit <- as.data.frame(
  reduce(l_results_params_gen, rbind)
  ) %>% as_tibble()
colnames(tbl_params_fit) <- c("c_fit", "w_fit", "bias_fit")

tbl_params_all <- cbind(
  tbl_params,
  tbl_params_fit
)


ggplot(tbl_params_all, aes(c, c_fit)) +
  geom_point()

tbl_params_all <- tbl_params_all %>%
  group_by(n_reps) %>%
  summarize(
    c_c = cor(c, c_fit),
    c_w = cor(c, w_fit),
    c_bias = cor(c, bias_fit),
    w_c = cor(w, c_fit),
    w_w = cor(w, w_fit),
    w_bias = cor(w, bias_fit),
    bias_c = cor(bias, c_fit),
    bias_w = cor(bias, w_fit),
    bias_bias = cor(bias, bias_fit)
  ) %>% pivot_longer(-n_reps) %>%
  mutate(
    gen = str_match(name, ("^(.+)_"))[, 2],
    recover = str_match(name, ("_(.+)$"))[, 2]
  )

ggplot(tbl_params_all, aes(gen, recover)) +
  geom_tile(aes(fill = value)) +
  geom_label(aes(label = str_c("r = ", round(value, 2)))) +
  facet_wrap(~ n_reps) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Generated", y = "Recovered") + 
  theme(strip.background = element_rect(fill = "white")) + 
  scale_fill_gradient2(low = "skyblue2", high = "tomato4", name = "Correlation")




# todo
# save tbl_samples_ii for recognition and RT analyses







