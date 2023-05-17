

library(tidyverse)
library(grid)
library(gridExtra)


# two category structures
# squared structure as in rep change project
# information integration condition using the identity function as decision boundary

# information integration
x1 <- seq(1, 10, by = 1)
x2 <- seq(1, 10, by = 1)
tbl_x_ii <- crossing(x1, x2)
tbl_x_ii$category <- tbl_x_ii$x1 < tbl_x_ii$x2
tbl_x_ii <- tbl_x_ii[!(tbl_x_ii$x1 == tbl_x_ii$x2), ]
tbl_x_ii$cat_structure <- "rule"
pl_rule <- ggplot(tbl_x_ii, aes(x1, x2, group = category)) + 
  geom_point(aes(color = category)) + 
  geom_abline() + 
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_viridis_d(name = "Category")
  

# squared category
tbl_x_sq <- crossing(x1, x2)
tbl_x_sq$category <- factor((tbl_x_sq$x1 > mean(x1)) + (tbl_x_sq$x2 > mean(x2)) * 2 )
tbl_x_sq$cat_structure <- "information-integration"

tbl_categories <- rbind(tbl_x_ii, tbl_x_sq)
pl_ii <- ggplot(tbl_x_sq, aes(x1, x2, group = category)) + 
  geom_point(aes(color = category)) + 
  geom_hline(yintercept = mean(x2)) +
  geom_vline(xintercept = mean(x1)) + 
  theme_bw() +
  scale_x_continuous(expand = c(0.1, 0.1), breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(expand = c(0.1, 0.1), breaks = seq(0, 10, by = 2)) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_viridis_d(name = "Category")

grid.draw(arrangeGrob(pl_rule, pl_ii, nrow = 1))






