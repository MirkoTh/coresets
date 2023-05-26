gcm_base <- function(x_new, tbl_x, n_feat, c, w, d_measure = 1){
  #' compute class probabilities with the GCM model
  #' 
  #' @description summed similarity computation with gcm;
  #' using sensitivity and attentional weighting;
  #' currently, 
  #' response bias is not implemented
  #' @param x_new the x coordinates of the new item
  #' @param tbl_x the tbl with all memory exemplars
  #' @param n_feat number of feature dimensions
  #' @param c sensitivity
  #' @param w attentional weighting
  #' @return a vector with the class probabilities for the new item
  l_x_cat <- split(tbl_x, tbl_x$category)
  sims_cat <- map(l_x_cat, f_similarity_cat, w, c, x_new, d_measure)
  map_dbl(sims_cat, ~ sum(.x)/sum(map_dbl(sims_cat, sum)))
  
}

f_similarity <- function(x1, x2, w, c, x_new, d_measure) {
  d <- ((w[1]*abs((x_new$x1 - x1))^d_measure + w[2]*abs((x_new$x2 - x2))^d_measure))^(1/d_measure)
  exp(-d*c)
}
f_similarity_cat <- function(x, w, c, x_new, d_measure) {
  pmap_dbl(x[, c("x1", "x2")], f_similarity, w, c, x_new, d_measure)
}



x_new <- tibble(x1 = 4.5, x2 = 4.7)
tbl_x <- tbl_x_ii_inb
n_feat <- 2
c <- .5
w <- rep(1/n_feat, n_feat) # equal
d_measure <- 1

gcm_base(x_new, tbl_x, 2, c, w, d_measure)