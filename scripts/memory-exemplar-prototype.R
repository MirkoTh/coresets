library(tidyverse)

v_time_passed <- seq(1, 10, by = 1)
v_repetition <- seq(1, 3, by = 1)
v_time_to_repetition <- c(3, 6, 9)

tbl_memory <- crossing(
  retention_interval = v_time_passed,
  n_repetition = v_repetition,
  time_spacing = v_time_to_repetition
)

tbl_memory <- tbl_memory %>%
  mutate(
    p_retention = ifelse(retention_interval ^ -2 == Inf, 1, retention_interval ^ -2),
    p_repetition = 1 - n_repetition ^ -.5,
    p_spacing = .075 * exp(-abs(time_spacing - retention_interval)*.8),
    p_memory = p_retention + p_repetition + p_spacing,
    p_memory = ifelse(p_memory > 1, 1, p_memory),
    weight_items = p_memory / (1 + p_memory), # assuming prototype is remembered perfectly
    time_spacing = as.factor(time_spacing),
    n_repetition = str_c("Nr. Repetitions = ", n_repetition)
  )

ggplot(tbl_memory, aes(retention_interval, p_memory, group = time_spacing)) +
  geom_line(aes(color = time_spacing)) +
  geom_point(color = "white", size = 3) +
  geom_point(aes(color = time_spacing)) + 
  theme_bw() +
  facet_wrap(~ n_repetition, nrow = 3) +
  scale_color_viridis_d(name = "Temporal Spacing") +
  scale_x_continuous(breaks = seq(2, 10, by = 2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Retention Interval", 
    y = "p (Retrieval)"
  ) + 
  theme_bw() +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme(strip.background = element_rect(fill = "white"))
  

ggplot(tbl_memory, aes(retention_interval, weight_items, group = time_spacing)) +
  geom_line(aes(color = time_spacing)) +
  geom_point(color = "white", size = 3) +
  geom_point(aes(color = time_spacing)) + 
  theme_bw() +
  facet_wrap(~ n_repetition, nrow = 1) +
  scale_color_viridis_d(name = "Temporal Spacing") +
  scale_x_continuous(breaks = seq(2, 10, by = 2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Retention Interval", 
    y = "Weight Items"
  )

