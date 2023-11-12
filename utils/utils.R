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
  x$lag <- abs(x$trial_id - max(x$trial_id)) + x_new$trial_id
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


category_probs <- function(x, tbl_transfer, tbl_x, n_feat, d_measure, lo, hi) {
  #' @description calculate category probabilities for every stimulus in the transfer set
  #' @param x parameters
  #' @param tbl_transfer transfer/test data
  #' @param tbl_x training data
  #' @param n_feat number of features
  #' @param d_measure distance measure, 1 for city-block, 2 for euclidean
  #' @param lo vector with lower bounds of parameters
  #' @param hi vector with upper bounds of parameters
  #' @return negative 2 * summed log likelihood
  #' 
  x <- pmap(list(x, lo, hi), upper_and_lower_bounds_revert)
  l_transfer_x <- split(tbl_transfer[, c("x1", "x2", "trial_id")], 1:nrow(tbl_transfer))
  l_category_probs <- map(
    l_transfer_x, gcm_base, 
    tbl_x = tbl_x, n_feat = n_feat, 
    c = x[["c"]], w = x[["w"]], bias = x[["bias"]], 
    delta = ifelse(is.null(x[["delta"]]), 0, x[["delta"]]), 
    d_measure = d_measure
  )
  tbl_probs <- as.data.frame(reduce(l_category_probs, rbind)) %>% mutate(response = tbl_transfer$response)
  tbl_probs$prob_correct <- pmap_dbl(
    tbl_probs[, c("0", "1", "response")],
    ~ c(..1, ..2)[as.numeric(as.character(..3)) + 1]
  )
  return(tbl_probs)
}


gcm_likelihood_no_forgetting <- function(x, tbl_transfer, tbl_x, n_feat, d_measure, lo, hi) {
  #' @description -2 * negative log likelihood of transfer set given training data
  #' and a gcm without a forgetting parameter (i.e., forgetting set to 0)
  #' @param x parameters
  #' @param tbl_transfer transfer/test data
  #' @param tbl_x training data
  #' @param n_feat number of features
  #' @param d_measure distance measure, 1 for city-block, 2 for euclidean
  #' @param lo vector with lower bounds of parameters
  #' @param hi vector with upper bounds of parameters
  #' @return negative 2 * summed log likelihood
  #' 
  tbl_probs <- category_probs(x, tbl_transfer, tbl_x, n_feat, d_measure, lo, hi)
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
  #' @param lo vector with lower bounds of parameters
  #' @param hi vector with upper bounds of parameters
  #' @return negative 2 * summed log likelihood
  #' 
  tbl_probs <- category_probs(x, tbl_transfer, tbl_x, n_feat, d_measure, lo, hi)
  ll <- log(tbl_probs$prob_correct)
  neg2llsum <- -2 * sum(ll)
  return(neg2llsum)
}


add_sample <- function(tbl_new_sample, tbl_base, params, tbl_transfer, n_feat, d_measure, f_likelihood, lo, hi) {
  #' @description evaluate likelihood on transfer/test set when a new sample is added to the training data
  f_likelihood(
    params, tbl_transfer, 
    tbl_x = rbind(tbl_base, tbl_new_sample) %>% mutate(trial_id = 1:(nrow(tbl_base) + 1)),
    n_feat = n_feat, d_measure = d_measure, lo = lo, hi = hi
  )
}


remove_sample <- function(rwn_remove, tbl_base, params, tbl_transfer, n_feat, d_measure, f_likelihood, lo, hi) {
  #' @description evaluate likelihood on transfer/test set when one training example is removed
  f_likelihood(
    params, tbl_transfer, tbl_x = tbl_base[-rwn_remove, ] %>% mutate(trial_id = 1:(nrow(tbl_base) - 1)),
    n_feat = n_feat, d_measure = d_measure, lo = lo, hi = hi
  )
}

remove_sample_general <- function(rwn_remove, tbl_base, m, tbl_transfer) {
  #' @description evaluate likelihood on transfer/test set when one training example is removed
  #' @param rwn_remove the row from the tbl_base to remove
  #' @param m the model including details about it
  #' @param tbl_transfer the data set to evaluate the model upon
  if (m$name == "gcm") {
    neg2loglik <- m$f_likelihood(
      m$params, tbl_transfer, tbl_x = tbl_base[-rwn_remove, ] %>% mutate(trial_id = 1:(nrow(tbl_base) - 1)),
      n_feat = m$n_feat, d_measure = m$d_measure, lo = m$lo, hi = m$hi
    )
  } else if (m$name == "svm") {
    # re-fit svm
    m$model <- svm(
      response ~ x1 + x2, data = tbl_base[-rwn_remove, ], type = "C-classification", 
      kernel = "linear", probability = TRUE
    )
    y_preds <- predict(m$model, tbl_transfer, probability = TRUE)
    lik <- pmap_dbl(cbind(attr(y_preds, "probabilities") %>% as.data.frame(), tbl_transfer$category), ~ c(..1, ..2)[..3])
    neg2loglik <- -2*sum(log(lik))
  }
  return(neg2loglik)
}



importance_upsampling <- function(tbl_importance, tbl_imb_plus, params, tbl_transfer, n_feat, d_measure, lo, hi, n_add = 10) {
  #' @description add n_add most important upsampled data points to the training set
  #' @param tbl_importance the data points to upsample from
  
  
  l_new_samples <- split(tbl_importance %>% dplyr::select(-trial_id), tbl_importance$trial_id)
  
  future::plan(multisession, workers = future::availableCores() - 2)
  l_samples <- l_new_samples
  l_tbl_add <- list()
  # sequential importance sampling
  # optimizes labels on training data
  for (i in 1:n_add) {
    v_importance <- future_map_dbl(
      l_samples, add_sample, tbl_base = tbl_imb_plus, params = params, 
      tbl_transfer = tbl_transfer, n_feat = 2, d_measure = 1,
      f_likelihood = gcm_likelihood_no_forgetting, lo = lo, hi = hi
    )
    tbl_importance$importance <- v_importance
    # pick best imagined data point
    tbl_importance <- tbl_importance %>% mutate(rank_importance = rank(importance, ties.method = "min"))
    tbl_imb_plus <- rbind(tbl_imb_plus, tbl_importance %>% dplyr::filter(rank_importance == 1) %>% dplyr::select(x1, x2, category, response))
    # exclude sampled points from sampling set
    idx_chosen_point <- which.min(v_importance)
    tbl_importance <- tbl_importance[-idx_chosen_point,]
    l_samples <- l_samples[-idx_chosen_point]
    l_tbl_add[[i]] <- tbl_imb_plus
  }
  
  future::plan("default")
  return(l_tbl_add)
}


importance_downsampling <- function(tbl_drop, m, tbl_transfer, cat_down, n_keep_max = 10) {
  #' @description remove n_keep least important upsampled data points from the training set
  
  future::plan(multisession, workers = future::availableCores() - 2)
  # sequential importance sampling
  
  l_tbl_drop <- list()
  
  n_start <- nrow(tbl_drop)
  n_keep_min <- 0
  nrow_drop <- nrow(tbl_drop %>% filter(category %in% cat_down))
  while (nrow_drop > n_keep_min) {
    rows_to_drop <- which(tbl_drop$category %in% cat_down)
    v_importance <- future_map_dbl(
      rows_to_drop, remove_sample_general, tbl_base = tbl_drop, m = m, 
      tbl_transfer = tbl_transfer,
      .options = furrr_options(seed = TRUE)
    )
    tbl_drop$importance <- max(v_importance + 1)
    tbl_drop$importance[rows_to_drop] <- v_importance
    # pick best imagined data point
    tbl_drop <- tbl_drop %>% mutate(rank_importance = rank(importance, ties.method = "min"))
    # we want to throw out that data point, which changes the -2LL the least
    # the lowest -2LL refers to the best prediction after exclusion of one data point
    idx_chosen_point <- which.min(tbl_drop$importance)
    tbl_drop <- tbl_drop[-idx_chosen_point,]
    nrow_drop <- nrow(tbl_drop %>% filter(category %in% cat_down))
    if (nrow_drop <= n_keep_max){
      l_tbl_drop[[nrow_drop]] <- tbl_drop
    }
    if (nrow_drop == 1) break
    cat(str_c("iteration ", n_start - nrow(tbl_drop), "\n"))
  }
  
  future::plan("default")
  return(l_tbl_drop)
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
    left_join(tbl_new[, c("x1", "x2", "is_new")], by = c("x1", "x2")) %>%
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


plot_new_and_dropped <- function(tbl_all) {
  ggplot(tbl_all, aes(x1, x2, group = category)) + 
    geom_point(aes(color = category, alpha = is_dropped)) + 
    geom_point(aes(size = as.numeric(is_new), alpha = is_dropped), shape = 1, color = "grey70") +
    geom_point(aes(size = as.numeric(is_dropped)), shape = 1, color = "grey70") +
    theme_bw() +
    scale_alpha_manual(values = c(1, 0), guide = "none") +
    scale_size_manual(guide = "none") +
    scale_x_continuous(expand = c(.1, .1), breaks = seq(2, 10, by = 2)) +
    scale_y_continuous(expand = c(.1, .1), breaks = seq(2, 10, by = 2)) +
    labs(x = expression(x[1]), y = expression(x[2])) +
    theme(strip.background = element_rect(fill = "white")) +
    scale_color_viridis_d(name = "Category") + 
    geom_abline() +
    scale_size_continuous(guide = "none") +
    ggtitle("Up- and Downsampling")
}


generate_data <- function(n_reps, params, tbl_transfer, tbl_train, n_feat, d_measure, lo, hi) {
  #' @description helper function to generate data on transfer data given 
  #' training data and set of parameters
  
  tbl_cat_probs <- category_probs(params$tf, tbl_transfer, tbl_train, 2, 1, lo, hi)
  tbl_generate <- tibble(
    repeat_tibble(tbl_transfer, n_reps), 
    prob_correct = rep(tbl_cat_probs$prob_correct, n_reps)
  )
  tbl_generate$accuracy <- rbernoulli(nrow(tbl_generate), tbl_generate$prob_correct)
  tbl_generate$response <- pmap_dbl(
    tbl_generate[, c("category", "accuracy")], 
    ~ c((as.numeric(as.character(..1)) - 1) * -1, as.numeric(as.character(..1)))[(..2 + 1)]
  )
  return(tbl_generate)
}


generate_and_fit <- function(n_reps, params, k, tbl_train_orig, l_tbl_train_strat, tbl_transfer, l_info, is_strategic) {
  #' @description generate data given model and strategically sampled data
  #' and fit them given strategically sampled data or originally presented data
  #' 
  #' n_feat, d_measure, lo, hi
  
  if (is_strategic) {
    # generate data given strat. sampling model
    tbl_train_strat <- l_tbl_train_strat[[k]] %>% mutate(trial_id = sample(1:nrow(.), nrow(.), replace = FALSE))
    tbl_generate <- generate_data(n_reps, params, tbl_transfer, tbl_train_strat, l_info$n_feat, l_info$d_measure, l_info$lo[1:3], l_info$hi[1:3])
  } else if (!is_strategic) {
    tbl_generate <- generate_data(n_reps, params, tbl_transfer, tbl_train_orig, l_info$n_feat, l_info$d_measure, l_info$lo, l_info$hi)
  }
  
  
  # starting values for strat. sampling and default gcm model
  params_init <- c(c = 1, w = .5, bias = .5)
  params_init_tf <- pmap(list(params_init[1:3], l_info$lo[1:3], l_info$hi[1:3]), upper_and_lower_bounds)
  # starting values for decay gcm model
  params_init_decay <- c(c = 1, w = .5, bias = .5, delta = 0.5)
  params_init_decay_tf <- pmap(list(params_init_decay, l_info$lo, l_info$hi), upper_and_lower_bounds)
  
  
  # iterate over all plausible ks
  l_params_strat <- list()
  l_results_strat <- list()
  for (i in 1:length(l_tbl_train_strat)) {
    results_strat <- optim(
      params_init_tf,
      gcm_likelihood_no_forgetting,
      tbl_transfer = tbl_generate,
      tbl_x = l_tbl_train_strat[[i]] %>% mutate(trial_id = sample(1:nrow(.), nrow(.), replace = FALSE)), 
      n_feat = l_info$n_feat,
      d_measure = l_info$d_measure,
      lo = l_info$lo[1:3],
      hi = l_info$hi[1:3]
    )
    
    params_strat <- list()
    params_strat[["not_tf"]] <- pmap_dbl(list(results_strat$par, lo[1:3], hi[1:3]), upper_and_lower_bounds_revert)
    params_strat[["tf"]] <- results_strat$par
    l_params_strat[[i]] <- params_strat
    l_results_strat[[i]] <- results_strat
  }
  
  
  results_orig <- optim(
    params_init_tf,
    gcm_likelihood_no_forgetting,
    tbl_transfer = tbl_generate,
    tbl_x = tbl_train_orig, 
    n_feat = l_info$n_feat,
    d_measure = l_info$d_measure,
    lo = l_info$lo[1:3],
    hi = l_info$hi[1:3]
  )
  
  params_orig <- list()
  params_orig[["not_tf"]] <- pmap_dbl(list(results_orig$par, l_info$lo[1:3], l_info$hi[1:3]), upper_and_lower_bounds_revert)
  params_orig[["tf"]] <- results_orig$par
  
  results_decay <- optim(
    params_init_decay_tf,
    gcm_likelihood_forgetting,
    tbl_transfer = tbl_generate,
    tbl_x = tbl_train_orig, 
    n_feat = l_info$n_feat,
    d_measure = l_info$d_measure,
    lo = l_info$lo,
    hi = l_info$hi
  )
  
  params_decay <- list()
  params_decay[["not_tf"]] <- pmap_dbl(list(results_decay$par, l_info$lo, l_info$hi), upper_and_lower_bounds_revert)
  params_decay[["tf"]] <- results_decay$par
  
  n2lls <- list()
  n2lls[["strategic"]] <- map(l_results_strat, "value")
  n2lls[["original"]] <- results_orig$value
  n2lls[["decay"]] <- results_decay$value
  
  return(list(
    n2lls = n2lls, 
    params_strat = l_params_strat, 
    params_orig = params_orig, 
    params_decay = params_decay
  ))
}


save_my_pdf_and_tiff <- function(pl, path_fl, w, h) {
  save_my_pdf(pl, str_c(path_fl, ".pdf"), w, h)
  save_my_tiff(pl, str_c(path_fl, ".tiff"), w, h)
}


sims_hotspot_strat <- function(w, sens, gamma, tbl_hotspots, tbl_strat, tbl_test) {
  #' @description calculate similarities of test stimuli to originally
  #' presented stimuli and to strategically sampled stimuli
  #' @param w attentional weighting parameter
  #' @param sens sensitivity parameter (aka c)
  #' @param gamma response scaling parameter
  #' @param tbl_hotspots presented data
  #' @param tbl_strat imagined data
  #' @param tbl_test test data
  #' @return a list with both similarities
  
  # similarity on presented data
  sims_hotspots <- pmap_dbl(
    tbl_test[, c("x1", "x2")], 
    ~ sum(pmap_dbl(
      tbl_hotspots[, c("x1", "x2")], 
      f_similarity, 
      c(w, 1 - w), sens, tibble(x1 = .x, x2 = .y), 1
    ))
  )
  
  # similarity on strategically sampled data
  sims_strat <- pmap_dbl(
    tbl_test[, c("x1", "x2")], 
    ~ sum(pmap_dbl(
      tbl_strat[, c("x1", "x2")], 
      f_similarity, 
      c(w, 1 - w), sens, tibble(x1 = .x, x2 = .y), 1
    ))
  )
  
  # response scaling and z scaling
  sims_strat_z <- scale(sims_strat ^ gamma)[, 1]
  sims_hotspots_z <- scale(sims_hotspots ^ gamma)[, 1]
  
  return(list(
    sims_strat_z = sims_strat_z, sims_hotspots_z = sims_hotspots_z
  ))
  
}


gen_wiener_2_slopes <- function(
    sim_hs, sim_strat, alpha, beta, tau, delta_ic, delta_sl1, delta_sl2, n_reps
) {
  #' @description generate recognition data with wiener likelihood 
  #' given similarities to presented and imagined data
  #' both similarities affect the drift rate according to 
  #' delta_sl1 and delta_sl2, respectively
  #' @return a tibble with generated data and similarities
  
  tbl_rt_gen <- as_tibble(map2_df(
    sim_hs, sim_strat, 
    ~ rwiener(
      n = n_reps, alpha = alpha, tau = tau, 
      beta = beta, delta = delta_ic + delta_sl1 *.x + delta_sl2 * .y
    )
  ))
  tbl_rt_gen$rt <- tbl_rt_gen$q
  tbl_rt_gen$model <- "Generating: 2 Slopes"
  tbl_rt_gen$resp_recode <- tbl_rt_gen$resp
  tbl_rt_gen$resp <- fct_relabel(tbl_rt_gen$resp, ~ c("old", "new"))
  tbl_rt_gen$sim_strat_z <- rep(sim_strat, each = n_reps)
  tbl_rt_gen$sim_hotspot_z <- rep(sim_hs, each = n_reps)
  
  return(tbl_rt_gen)
}

wiener_reg1_delta_log <- function(x, my_tbl) {
  #' @description Wiener LL with linear regression on drift rate
  #' only regressing drift rate on one similarity
  #' @return the -2 * sum of the LL
  
  
  alpha <- x[["alpha"]]
  tau <- x[["tau"]]
  beta <- x[["beta"]]
  delta_ic <- x[["delta_ic"]]
  delta_slope <- x[["delta_sl1"]]
  
  lik <- pmap_dbl(
    my_tbl[, c("rt", "resp_recode", "pred_lr")], ~ dwiener(
      q = ..1, alpha = alpha, tau = tau, beta = beta, 
      delta = delta_ic + delta_slope * ..3, resp = ..2
    )
  )
  
  neg2loglik <- -2*sum(log(pmax(lik,1e-10)))
  
  return(neg2loglik)
}

wiener_reg2_delta_log <- function(x, my_tbl) {
  #' @description Wiener LL with linear regression on drift rate
  #' regressing drift rate on two similarities
  #' @return the -2 * sum of the LL
  
  alpha <- x[["alpha"]]
  tau <- x[["tau"]]
  beta <- x[["beta"]]
  delta_ic <- x[["delta_ic"]]
  delta_slope1 <- x[["delta_sl1"]]
  delta_slope2 <- x[["delta_sl2"]]
  
  
  lik <- pmap_dbl(
    my_tbl[, c("rt", "resp_recode", "pred_lr1", "pred_lr2")], 
    ~ dwiener(
      q = ..1, alpha = alpha, tau = tau, beta = beta, 
      delta = delta_ic + delta_slope1 * ..3 + delta_slope2 * ..4,
      resp = ..2
    )
  )
  
  neg2loglik <- -2*sum(log(pmax(lik,1e-10)))
  
  return(neg2loglik)
}

gen_recognition_2_slopes <- function(params, tbl_hotspots, tbl_strat, tbl_test) {
  #' @description convenience function generating recognition responses
  #' given global similarity computation of presented and imagined data
  #' and wiener diffusion process while regressing drift rate on 
  #' these two similarities
  #' @return a tibble with the generated responses
  
  l_sims <- sims_hotspot_strat(
    params$w, params$sens, params$gamma, tbl_hotspots, tbl_strat, tbl_test
  )
  tbl_gen <- gen_wiener_2_slopes(
    l_sims$sims_hotspots_z, l_sims$sims_strat_z,
    params$alpha, params$beta, params$tau,
    params$delta_ic, params$delta_sl1, params$delta_sl2,
    1
  )
  tbl_gen$x1 <- tbl_test$x1
  tbl_gen$x2 <- tbl_test$x2
  return(tbl_gen)
}


ll_recognition_2_slopes <- function(x, tbl_hotspots, tbl_strat, tbl_test) {
  #' @description helper function to calculate neg 2 * LL fitting 
  #' parameters of representation model and decision model
  #' assuming two similarities affect drift rate
  #' @return the value of the -2 * LL
  #'
  x <- pmap(list(x, lo2, hi2), upper_and_lower_bounds_revert)
  l_sims <- sims_hotspot_strat(
    x[["w"]], x[["sens"]], x[["gamma"]], tbl_hotspots, tbl_strat, tbl_test
  )
  tbl_test$pred_lr1 <- l_sims$sims_hotspots_z
  tbl_test$pred_lr2 <- l_sims$sims_strat_z
  neg2ll <- wiener_reg2_delta_log(x, tbl_test)
  
  return(neg2ll)
}

ll_recognition_1_slope <- function(x, tbl_hotspots, tbl_strat, tbl_test, sim_which) {
  #' @description helper function to calculate neg 2 * LL fitting 
  #' parameters of representation model and decision model
  #' assuming only one similarity affects drift rate
  #' @return the value of the -2 * LL
  #'   
  x <- pmap(list(x, lo1, hi1), upper_and_lower_bounds_revert)
  l_sims <- sims_hotspot_strat(
    x[["w"]], x[["sens"]], x[["gamma"]], tbl_hotspots, tbl_strat, tbl_test
  )
  tbl_test$pred_lr <- l_sims[[sim_which]]
  neg2ll <- wiener_reg1_delta_log(x, tbl_test)
  
  return(neg2ll)
}


sample_from_grid <- function(tbl_x, n_trials_total, my_sd) {
  #' @description sample x from bivariate normal for n_trials_total times
  #' @return all sampled data points as a tbl
  #'  
  n_stim <- tbl_x %>% count(category)
  n_stimulus_reps <- tibble(
    category = n_stim$category,
    n = (n_trials_total / nrow(n_stim)) / n_stim$n
  )
  tbl_x <- tbl_x %>% left_join(n_stimulus_reps, by = "category")
  
  overlap <- 1
  # do not allow samples to come from wrong category
  while(overlap > 0) {
    l_samples <- pmap(tbl_x[, c("x1", "x2", "category", "n")], sample_2d_stimuli, sd = my_sd)
    tbl_samples <- reduce(l_samples, rbind)
    tbl_overlap <- tbl_samples %>% group_by(category) %>%
      count(cat_wrong_1 = x1 >= x2, cat_wrong_0 = x2 >= x1) %>%
      pivot_longer(cols = c(cat_wrong_0, cat_wrong_1)) %>%
      mutate(cat_wrong = str_extract(name, "[0-9]$")) %>%
      filter(value & as.numeric(as.character(category)) == cat_wrong)
    if (nrow(tbl_overlap) == 0) {overlap <- 0}
  }
  tbl_samples <- tbl_samples %>%
    mutate(trial_id = sample(1:nrow(.), nrow(.)))
}
