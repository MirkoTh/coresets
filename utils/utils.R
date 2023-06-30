gcm_base <- function(x_new, tbl_x, n_feat, c, w, bias, delta, d_measure = 1){
  #' compute class probabilities with the GCM model
  #' 
  #' @description summed similarity computation with gcm;
  #' using sensitivity, attentional weighting, and response bias;
  #' @param x_new the x coordinates of the new item
  #' @param tbl_x the tbl with all memory exemplars,
  #' including a column "category" denoting the category of that item
  #' @param n_feat number of feature dimensions
  #' @param c sensitivity
  #' @param w attentional weighting
  #' @param bias response bias / category prior
  #' @param delta forgetting rate (if delta == 0, no forgetting)
  #' @return a vector with the class probabilities for the new item
  w <- c(w, 1 - w)
  bias <- c(bias, 1 - bias)
  l_x_cat <- split(tbl_x, tbl_x$category)
  sims_cat <- map(l_x_cat, f_similarity_cat, w, c, delta, x_new, d_measure)
  sims_cat_sum <- map_dbl(sims_cat, sum)
  sims_cat_sum_biased <- sims_cat_sum * bias
  map_dbl(sims_cat_sum_biased, ~ .x/sum(sims_cat_sum_biased))
}

f_similarity <- function(x1, x2, w, c, x_new, d_measure) {
  #' @description helper function calculating similarity for one item
  #' given a currently presented item to be classified
  d <- (
    w[1]*abs((x_new$x1 - x1))^d_measure + w[2]*abs((x_new$x2 - x2))^d_measure
  )^(1/d_measure)
  exp(-d*c)
}

f_similarity_cat <- function(x, w, c, delta, x_new, d_measure) {
  #' @description helper function calculating similarities for all items
  #' within a given category
  x$lag <- abs(x$trial_id - max(x$trial_id))
  x$prop_decay <- exp(- (delta * x$lag))
  sims <- pmap_dbl(x[, c("x1", "x2")], f_similarity, w, c, x_new, d_measure)
  return(sims*x$prop_decay)
}

my_weighted_sample <- function(prop_1, tbl_df) {
  tibble(
    x1 = tbl_df$x1_1 * prop_1 + tbl_df$x1_2 * (1 - prop_1),
    x2 = tbl_df$x2_1 * prop_1 + tbl_df$x2_2 * (1 - prop_1)
  )
}

upsample <- function(category_used, n_upsample, tbl_x) {
  #' upsample by linearly interpolating between observed data points
  #' 
  #' @description divides up n_upsample equally by number of connections
  #' and linearly interpolates values with equal spacing
  #' @param category_used category to upsample data from
  #' @param n_upsample total number of data points to upsample
  #' @param tbl_x tbl_df with physically presented stimuli
  #' @return a tbl with the original data points enriched
  #' with the upsampled data points
  
  data_category <- tbl_x %>% dplyr::filter(category == category_used)
  n_unique_obs <- nrow(data_category)
  
  n_pairs <- ((n_unique_obs^2)/2) - (n_unique_obs/2)
  n_new_per_connection <- floor(n_upsample / n_pairs)
  n_new_per_connection <- ifelse(is.na(n_new_per_connection), 0, n_new_per_connection)
  
  v_range <- seq(0, 1, length.out = (n_new_per_connection + 2))
  v_range <- v_range[!(v_range %in% c(0, 1))]
  fully_crossed <- crossing(
    data_category[, c("x1", "x2")] %>% rename(x1_1 = x1, x2_1 = x2), 
    data_category[, c("x1", "x2")] %>% rename(x1_2 = x1, x2_2 = x2)
  )
  mat_help <- matrix(seq(1, nrow(fully_crossed), by = 1), nrow = n_unique_obs)
  unique_pairs <- fully_crossed[mat_help[lower.tri(mat_help)], ]
  
  unique_points <- data_category %>% dplyr::select(x1, x2)
  
  l_upsample <- map(v_range, my_weighted_sample, tbl_df = unique_pairs)
  if (length(l_upsample) > 0) {
    new_points <- reduce(l_upsample, rbind)
  } else if (length(l_upsample) == 0) {
    new_points <- tibble(x1 = numeric(), x2 = numeric())
  }
  tbl_upsample <- rbind(unique_points, new_points)
  tbl_upsample$category <- category_used
  
  return(tbl_upsample)
}





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



# todos
# add/remove an item from the set of presented items
# and test how much that changes the fit
# this makes sense only in the non-decay model
# because recently added items have more weight in the decay model


# number of finally used items is controlled by capacity parameter

tbl_x_ii_inb$response <- tbl_x_ii_inb$category

# category 0 example
gcm_base(tbl_x_ii_inb[1, c("x1", "x2")], tbl_x_ii_inb[-1, ], 2, c = .5, w = c(.5, .5), .5, 0, d_measure = 1)
# category 1 example
gcm_base(tbl_x_ii_inb[17, c("x1", "x2")], tbl_x_ii_inb, 2, c = 3, w = c(.001, .999), .5, 0, d_measure = 1)


# example using inbalanced data set sequentially presenting data

tbl_inb <- tbl_x_ii %>% 
  filter(
    category %in% c(0) |
      (category %in% c(1) & x1 %in% seq(1, 9, by = 4) & x2 %in% seq(2, 10, by = 4))
  )
plot_grid(tbl_inb) + geom_abline()

# for one category
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



# create 10 x 10 grid of data points
xs <- seq(1, 10, by = 1)
tbl_transfer <- crossing(x1 = xs, x2 = xs) %>% mutate(category = fct_rev(factor(x1 > x2, labels = c(1, 0))))
l_transfer_x <- split(tbl_transfer[, c("x1", "x2")], 1:nrow(tbl_transfer))
plot_grid(tbl_transfer) + geom_abline()

gcm_base(tbl_transfer[10, c("x1", "x2")], tbl_inb, n_feat = 2, c = 1, w = .5, bias = .5, delta = 0, d_measure = 1)


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
  scale_fill_gradient(name = "Prob. Error", low = "white", high = "tomato", limits = c(0, 1)) + 
  geom_abline() + ggtitle("Baseline Classification, Bias = .1")

# plot together
grid.draw(arrangeGrob(pl_inb, pl_gcm_baseline.5, pl_gcm_baseline.1, nrow = 1))



gcm_likelihood_no_forgetting <- function(x, tbl_transfer, tbl_x, n_feat, d_measure) {
  l_transfer_x <- split(tbl_transfer[, c("x1", "x2")], 1:nrow(tbl_transfer))
  l_category_probs <- map(l_transfer_x, gcm_base, tbl_x = tbl_x, n_feat = n_feat, c = x[[1]], w = x[[2]], bias = x[[3]], delta = 0, d_measure = d_measure)
  tbl_probs <- as.data.frame(reduce(l_category_probs, rbind)) %>% mutate(category = tbl_transfer$category)
  tbl_probs$prob_correct <- pmap_dbl(
    tbl_probs[, c("0", "1", "category")],
    ~ c(..1, ..2)[as.numeric(as.character(..3)) + 1]
  )
  ll <- log(tbl_probs$prob_correct)
  neg2llsum <- -2 * sum(ll)
  return(neg2llsum)
}

gcm_likelihood_forgetting <- function(x, tbl_transfer, tbl_x, n_feat, d_measure) {
  l_transfer_x <- split(tbl_transfer[, c("x1", "x2")], 1:nrow(tbl_transfer))
  l_category_probs <- map(l_transfer_x, gcm_base, tbl_x = tbl_x, n_feat = n_feat, c = x[[1]], w = x[[2]], bias = x[[3]], delta = x[[4]], d_measure = d_measure)
  tbl_probs <- as.data.frame(reduce(l_category_probs, rbind)) %>% mutate(category = tbl_transfer$category)
  tbl_probs$prob_correct <- pmap_dbl(
    tbl_probs[, c("0", "1", "category")],
    ~ c(..1, ..2)[as.numeric(as.character(..3)) + 1]
  )
  ll <- log(tbl_probs$prob_correct)
  neg2llsum <- -2 * sum(ll)
  return(neg2llsum)
}


# for one category
categories <- c(0, 1)
ns_upsample <- c(0, 100)

l_tbl_upsample_inb <- map2(categories, ns_upsample, upsample, tbl_x = tbl_inb)
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


add_sample <- function(tbl_new_sample, tbl_base, params, tbl_transfer, n_feat, d_measure) {
  gcm_likelihood_no_forgetting(params, tbl_transfer, tbl_x = rbind(tbl_base, tbl_new_sample) %>% mutate(trial_id = 1:(nrow(tbl_base) + 1)), n_feat = 2, d_measure = 1)
}




# varying number of upsampled points and evaluating likelihood of given responses
# possible solution: first, upsample sufficient number of data points (e.g., 50)
# then, do a constrained optimization with n ranging between 0 and that number (i.e., 50)

importance_upsampling <- function(l_new_samples, tbl_importance, tbl_inb_plus, params, tbl_transfer, n_feat, d_measure, n_max = 10) {
  future::plan(multisession, workers = future::availableCores() - 2)
  l_samples <- l_new_samples
  # sequential importance sampling
  for (i in 1:n_max) {
    v_importance <- future_map_dbl(l_samples, add_sample, tbl_base = tbl_inb_plus, params = params, tbl_transfer = tbl_transfer, n_feat = 2, d_measure = 1)
    tbl_importance$importance <- v_importance
    # pick best imagined data point
    tbl_importance <- tbl_importance %>% mutate(rank_importance = rank(importance, ties.method = "min"))
    tbl_inb_plus <- rbind(tbl_inb_plus, tbl_importance %>% dplyr::filter(rank_importance == 1) %>% dplyr::select(x1, x2, category))
    # exclude sampled points from sampling set
    idx_chosen_point <- which.min(v_importance)
    tbl_importance <- tbl_importance[-idx_chosen_point,]
    l_samples <- l_samples[-idx_chosen_point]
  }
  
  future::plan("default")
  return(tbl_inb_plus)
}
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


remove_sample <- function(rwn_remove, tbl_base, params, tbl_transfer, n_feat, d_measure) {
  gcm_likelihood_no_forgetting(params, tbl_transfer, tbl_x = tbl_base[-rwn_remove, ] %>% mutate(trial_id = 1:(nrow(tbl_base) - 1)), n_feat = 2, d_measure = 1)
}



importance_downsampling <- function(tbl_inb_plus, params, tbl_transfer, n_feat, d_measure, n_max = 10) {
  future::plan(multisession, workers = future::availableCores() - 2)
  # sequential importance sampling
  tbl_drop <- tbl_inb_plus
  while (nrow(tbl_drop) > n_max) {
    v_importance <- future_map_dbl(1:nrow(tbl_drop), remove_sample, tbl_base = tbl_drop, params = params, tbl_transfer = tbl_transfer, n_feat = 2, d_measure = 1)
    tbl_drop$importance <- v_importance
    # pick best imagined data point
    tbl_drop <- tbl_drop %>% mutate(rank_importance = rank(importance, ties.method = "min"))
    # exclude sampled points from sampling set
    idx_chosen_point <- which.min(v_importance)
    tbl_drop <- tbl_drop[-idx_chosen_point,]
  }
  
  future::plan("default")
  return(tbl_drop)
}

plot_grid(tbl_drop)




