rm(list = ls())

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

