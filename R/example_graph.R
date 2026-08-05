library(tidyverse)
library(latex2exp)
source("R/main_theme.R")


example_graph <- function(S1, S2, a1, a2, b1, b2) {

tex_label <- TeX(sprintf(r"($f(t)=\frac{%s}{1+e^{%s(%s - t)}} + \frac{%s}{1+e^{%s(%s - t)}}$)", S1, a1, b1, S2, a2, b2), output = "character")

df_example <- tibble(time = seq(0, 48, 0.25)) %>%
  mutate(
    growth = S1 / (1 + exp(a1 * (b1 - time))),
    decay = S2 / (1 + exp(a2 * (b2 - time))),
    combined = growth + decay,
    raw = combined + rnorm(length(time), 0, 0.6)
  )

df_example %>%
  pivot_longer(c(growth, decay, combined), names_to = "type", values_to = "value") %>%
  mutate(type = factor(type, levels = c("growth", "decay", "combined"), labels = c("Primary", "Secondary", "Combined"))) %>%
  ggplot(aes(time, value, color = type)) +
  geom_point(aes(y = raw), size = 0.5, color = "#adafae") +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = color_palette) +
  annotate(
    "text",
    x = 20, y = 5, hjust = 0, size = 6, parse = TRUE,
    label = tex_label, color = "#adafae"
  ) +
  labs(x = "Time (hr)", y = "Normalized Fluorescence") +
  dark_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0, .95),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
  )

}