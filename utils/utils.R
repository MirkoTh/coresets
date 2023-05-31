gcm_base <- function(x_new, tbl_x, n_feat, c, w, delta, d_measure = 1){
  #' compute class probabilities with the GCM model
  #' 
  #' @description summed similarity computation with gcm;
  #' using sensitivity and attentional weighting;
  #' currently, 
  #' response bias is not implemented
  #' @param x_new the x coordinates of the new item
  #' @param tbl_x the tbl with all memory exemplars,
  #' including a column "category" denoting the category of that item
  #' @param n_feat number of feature dimensions
  #' @param w attentional weighting
  #' @param c sensitivity
  #' @param delta forgetting rate (if delta == 0, no forgetting)
  #' @return a vector with the class probabilities for the new item
  l_x_cat <- split(tbl_x, tbl_x$category)
  sims_cat <- map(l_x_cat, f_similarity_cat, w, c, delta, x_new, d_measure)
  map_dbl(sims_cat, ~ sum(.x)/sum(map_dbl(sims_cat, sum)))
  
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

upsample <- function(tbl_x, category_used, n_upsample) {
  #' upsample by linearly interpolating between observed data points
  #' 
  #' @description divides up n_upsample equally by number of connections
  #' and linearly interpolates values with equal spacing
  #' @param tbl_x tbl_df with physically presented stimuli
  #' @param category_used category to upsample data from
  #' @param n_upsample total number of data points to upsample
  #' @return a tbl with the original data points enriched
  #' with the upsampled data points
  
  data_category <- tbl_x %>% dplyr::filter(category == category_used)
  n_unique_obs <- nrow(data_category)
  
  n_pairs <- ((n_unique_obs^2)/2) - (n_unique_obs/2)
  n_new_per_connection <- floor(n_upsample / n_pairs)
  n_new_per_connection <- ifelse(is.na(n_new_per_connection), 0, n_new_per_connection)
  
  v_range <- seq(0, 1, length.out = (n_new_per_connection + 2))
  v_range <- v_range[2 : (length(v_range) - 1)]
  fully_crossed <- crossing(
    data_category[, c("x1", "x2")] %>% rename(x1_1 = x1, x2_1 = x2), 
    data_category[, c("x1", "x2")] %>% rename(x1_2 = x1, x2_2 = x2)
  )
  mat_help <- matrix(seq(1, nrow(fully_crossed), by = 1), nrow = n_unique_obs)
  unique_pairs <- fully_crossed[mat_help[lower.tri(mat_help)], ]

  unique_points <- data_category %>% dplyr::select(x1, x2)
  
  l_upsample <- map(v_range, my_weighted_sample, tbl_df = unique_pairs)
  tbl_upsample <- rbind(unique_points, reduce(l_upsample, rbind))
  tbl_upsample$category <- category_used
  # tbl_upsample <- tbl_upsample %>%
  #   group_by(x1, x2, category) %>%
  #   count() %>% dplyr::select(-n)
  
  return(tbl_upsample)
}





x_new <- tibble(x1 = 4.5, x2 = 4.7)
tbl_x <- tbl_x_ii_inb
n_feat <- 2
c <- .5
w <- rep(1/n_feat, n_feat) # equal
delta <- .99
d_measure <- 1

gcm_base(x_new, tbl_x, 2, c, w, 0, d_measure)


plot_grid(tbl_x) +
  geom_abline() + coord_cartesian(xlim = c(1, 10), ylim = c(1, 10))

# for one category
category_used <- 0
n_items_in_category <- nrow(tbl_x[tbl_x$category == category_used, ])
n_upsample <- 1000
tbl_upsample <- upsample(tbl_x, category_used, n_upsample)
plot_grid(tbl_upsample %>% mutate(category = factor(category)))

gcm_base(x_new, tbl_x, 2, c, w, delta, d_measure)


# todos
# add/remove an item from the set of presented items
# and test how much that changes the fit
# this makes sense only in the non-decay model
# because recently added items have more weight in the decay model


# number of finally used items is controlled by capacity parameter


















