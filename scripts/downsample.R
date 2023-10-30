# Load packages and utils -------------------------------------------------


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

# training set
# transfer set = training set
# upsampling set
# set of parameters, possibly fitted


# generate information integration data set
x1 <- seq(1, 6, by = 1)
x2 <- seq(1, 6, by = 1)
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
n_trials_total <- 300
tbl_samples_ii <- sample_from_grid(tbl_x_ii, n_trials_total)

plot_grid(tbl_samples_ii) + geom_abline()

tbl_transfer <- tbl_samples_ii
tbl_samples_ii <- simulate_responses(tbl_samples_ii)
tbl_transfer <- simulate_responses(tbl_transfer)
















