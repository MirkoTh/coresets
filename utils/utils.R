plot_grid <- function(tbl) {
  #' @description plot 2D grid of x1 and x2 data points
    ggplot(tbl, aes(x1, x2, group = category)) + 
    geom_point(aes(color = category)) + 
    theme_bw() +
    scale_x_continuous(expand = c(.1, .1), breaks = seq(2, 10, by = 2)) +
    scale_y_continuous(expand = c(.1, .1), breaks = seq(2, 10, by = 2)) +
    labs(x = expression(x[1]), y = expression(x[2])) +
    theme(strip.background = element_rect(fill = "white")) +
    scale_color_viridis_d(name = "Category")
}



sample_2d_stimuli <- function(x1, x2, category, n, sd) {
  #' @description sample from bivariate normal
  #' @param x1 mean 1
  #' @param x2 mean 2
  #' @param category a category label
  #' @param n number of samples
  #' @param sd the sd
  s <- mvrnorm(n, c(x1, x2), matrix(c(sd, 0, 0, sd), nrow = 2))
  tbl <- tibble(as.data.frame(s), category)
  colnames(tbl) <- c("x1", "x2", "category")
  return(tbl)
}


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
  #' @param d_measure distance measure, 1 for city-block, 2 for euclidean
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


gcm_likelihood_no_forgetting <- function(x, tbl_transfer, tbl_x, n_feat, d_measure, lo, hi) {
  #' @description -2 * negative log likelihood of transfer set given training data
  #' and a gcm without a forgetting parameter (i.e., forgetting set to 0)
  #' @param x parameters
  #' @param tbl_transfer transfer/test data
  #' @param tbl_x training data
  #' @param n_feat number of features
  #' @param d_measure distance measure, 1 for city-block, 2 for euclidean
  #' @return negative 2 * summed log likelihood
  #' 
  x <- pmap(list(x, lo, hi), upper_and_lower_bounds_revert)
  l_transfer_x <- split(tbl_transfer[, c("x1", "x2")], 1:nrow(tbl_transfer))
  l_category_probs <- map(l_transfer_x, gcm_base, tbl_x = tbl_x, n_feat = n_feat, c = x[[1]], w = x[[2]], bias = x[[3]], delta = 0, d_measure = d_measure)
  tbl_probs <- as.data.frame(reduce(l_category_probs, rbind)) %>% mutate(response = tbl_transfer$response)
  tbl_probs$prob_correct <- pmap_dbl(
    tbl_probs[, c("0", "1", "response")],
    ~ c(..1, ..2)[as.numeric(as.character(..3)) + 1]
  )
  ll <- log(tbl_probs$prob_correct)
  neg2llsum <- -2 * sum(ll)
  return(neg2llsum)
}


gcm_likelihood_forgetting <- function(x, tbl_transfer, tbl_x, n_feat, d_measure, lo, hi) {
  #' @description -2 * negative log likelihood of transfer set given training data
  #' and a gcm with a forgetting parameter
  #' @param x parameters
  #' @param tbl_transfer transfer/test data
  #' @param tbl_x training data
  #' @param n_feat number of features
  #' @param d_measure distance measure, 1 for city-block, 2 for euclidean
  #' @return negative 2 * summed log likelihood
  #' 
  x <- pmap(list(x, lo, hi), upper_and_lower_bounds_revert)
  l_transfer_x <- split(tbl_transfer[, c("x1", "x2")], 1:nrow(tbl_transfer))
  l_category_probs <- map(l_transfer_x, gcm_base, tbl_x = tbl_x, n_feat = n_feat, c = x[[1]], w = x[[2]], bias = x[[3]], delta = x[[4]], d_measure = d_measure)
  tbl_probs <- as.data.frame(reduce(l_category_probs, rbind)) %>% mutate(response = tbl_transfer$response)
  tbl_probs$prob_correct <- pmap_dbl(
    tbl_probs[, c("0", "1", "response")],
    ~ c(..1, ..2)[as.numeric(as.character(..3)) + 1]
  )
  ll <- log(tbl_probs$prob_correct)
  neg2llsum <- -2 * sum(ll)
  return(neg2llsum)
}


add_sample <- function(tbl_new_sample, tbl_base, params, tbl_transfer, n_feat, d_measure, f_likelihood, lo, hi) {
  #' @description evaluate likelihood on transfer/test set when a new sample is added to the training data
  f_likelihood(
    params, tbl_transfer, 
    tbl_x = rbind(tbl_base, tbl_new_sample) %>% mutate(trial_id = 1:(nrow(tbl_base) + 1)),
    n_feat = 2, d_measure = 1, lo = lo, hi = hi
    )
}


remove_sample <- function(rwn_remove, tbl_base, params, tbl_transfer, n_feat, d_measure, f_likelihood, lo, hi) {
  #' @description evaluate likelihood on transfer/test set when one training example is removed
  f_likelihood(
    params, tbl_transfer, tbl_x = tbl_base[-rwn_remove, ] %>% mutate(trial_id = 1:(nrow(tbl_base) - 1)),
    n_feat = 2, d_measure = 1, lo = lo, hi = hi
    )
}


importance_upsampling <- function(tbl_importance, tbl_inb_plus, params, tbl_transfer, n_feat, d_measure, lo, hi, n_max = 10) {
  #' @description add n_max most important upsampled data points to the training set

  
  l_new_samples <- split(tbl_importance %>% dplyr::select(-trial_id), tbl_importance$trial_id)
  
  future::plan(multisession, workers = future::availableCores() - 2)
  l_samples <- l_new_samples
  # sequential importance sampling
  # optimizes labels on training data
  for (i in 1:n_max) {
    v_importance <- future_map_dbl(
      l_samples, add_sample, tbl_base = tbl_inb_plus, params = params, 
      tbl_transfer = tbl_transfer, n_feat = 2, d_measure = 1,
      f_likelihood = gcm_likelihood_no_forgetting, lo = lo, hi = hi
      )
    tbl_importance$importance <- v_importance
    # pick best imagined data point
    tbl_importance <- tbl_importance %>% mutate(rank_importance = rank(importance, ties.method = "min"))
    tbl_inb_plus <- rbind(tbl_inb_plus, tbl_importance %>% dplyr::filter(rank_importance == 1) %>% dplyr::select(x1, x2, category, response))
    # exclude sampled points from sampling set
    idx_chosen_point <- which.min(v_importance)
    tbl_importance <- tbl_importance[-idx_chosen_point,]
    l_samples <- l_samples[-idx_chosen_point]
  }
  
  future::plan("default")
  return(tbl_inb_plus)
}


importance_downsampling <- function(tbl_inb, params, tbl_transfer, n_feat, d_measure, lo, hi, cat_down, n_max = 10) {
  #' @description remove n_max least important upsampled data points from the training set

  future::plan(multisession, workers = future::availableCores() - 2)
  # sequential importance sampling
  tbl_drop <- tbl_inb
  
  
  while (nrow(tbl_drop %>% filter(category == cat_down)) > n_max) {
    rows_to_drop <- which(tbl_drop$category == cat_down)
    v_importance <- future_map_dbl(
      rows_to_drop, remove_sample, tbl_base = tbl_drop, params = params, 
      tbl_transfer = tbl_transfer, n_feat = 2, d_measure = 1,
      f_likelihood = gcm_likelihood_no_forgetting, lo = lo, hi = hi
      )
    tbl_drop$importance <- max(v_importance + 1)
    tbl_drop$importance[rows_to_drop] <- v_importance
    # pick best imagined data point
    tbl_drop <- tbl_drop %>% mutate(rank_importance = rank(importance, ties.method = "min"))
    # exclude sampled points from sampling set
    idx_chosen_point <- which.min(tbl_drop$importance)
    tbl_drop <- tbl_drop[-idx_chosen_point,]
  }
  
  future::plan("default")
  return(tbl_drop)
}


upper_and_lower_bounds <- function(par, lo, hi) {
  log(((par - lo) / (hi - lo)) / (1 - (par - lo) / (hi - lo)))
}

upper_and_lower_bounds_revert <- function(par, lo, hi) {
  lo + ((hi - lo) / (1 + exp(-par)))
}



repeat_tibble <- function(tbl_df, n_reps) {
  #' concatenate the same tibble several times
  #' 
  #' @description copy a tibble n_reps times and rbind it to the original tibble
  #' @param tbl_df the tbl to be repeated
  #' @param n_reps the number of times the tbl should be repeated
  #' @return the new larger tibble
  
  i <- 1
  tbl_df_new <- tbl_df
  while (i < n_reps) {
    tbl_df_new <- rbind(tbl_df_new, tbl_df)
    i <- i + 1
  }
  return(tbl_df_new)
}



simulate_responses <- function(tbl_df) {
  #' @description simulate responses for two categories with diagonal 
  #' information integration structure
  #' @param tbl_df the tbl with x1 and x2 coordinates
  #' @return the same tbl with p_correct, accuracy, and responses added as columns


  tbl_df$x1_bd <- (tbl_df$x1 + tbl_df$x2) / 2
  tbl_df$x2_bd <- tbl_df$x1_bd
  tbl_df$d_bd <- sqrt((tbl_df$x1 - tbl_df$x1_bd) ^ 2 + (tbl_df$x2 - tbl_df$x2_bd) ^2)
  tbl_df$p_correct <- 1 - exp(-tbl_df$d_bd * 1.5)
  tbl_df$accuracy <- rbernoulli(nrow(tbl_df), p = tbl_df$p_correct)
  tbl_df$response <- pmap_dbl(
    tbl_df[, c("category", "accuracy")], 
    ~ ifelse(..2, as.numeric(as.character(..1)), abs(as.numeric(as.character(..1)) - 1))
  )
  return(tbl_df)
}



add_jitter <- function(tbl_df) {
  #' @description simulate responses for two categories with diagonal 
  #' information integration structure
  #' @param tbl_df the tbl with x1 and x2 coordinates
  #' @return the same tbl with p_correct, accuracy, and responses added as columns

  tbl_df %>% dplyr::select(x1, x2, category) %>%
    mutate(
      x1 = x1 + rnorm(nrow(.), 0, .01),
      x2 = x2 + rnorm(nrow(.), 0, .01)
    )
}


mark_changes <- function(tbl_final, tbl_original) {
  #' @description mark up- and downsampled points
  #' @param tbl_final the remaining tbl after strategic sampling
  #' @param tbl_original the tbl with all points presented to participants
  #' @return a tbl including all up- and downsampled points
  
  tbl_new <- tbl_final %>% 
    left_join(tbl_original[, c("x1", "x2", "cat_structure")], by = c("x1", "x2")) %>%
    mutate(is_new = is.na(cat_structure)) %>%
    dplyr::select(-cat_structure)
  # mark downsampled points
  tbl_dropped <- tbl_original %>% 
    left_join(tbl_final[, c("x1", "x2", "is_new")], by = c("x1", "x2")) %>%
    mutate(is_dropped = is.na(is_new)) %>%
    dplyr::select(-is_new)
  # integrate new and dropped points
  full_join(
    tbl_new[, c("x1", "x2", "category", "is_new")], 
    tbl_dropped[, c("x1", "x2", "category", "is_dropped")], 
    by = c("x1", "x2", "category")
  ) %>%
    replace_na(list(is_new = FALSE, is_dropped = FALSE))
}