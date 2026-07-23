main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 20),
  legend.title = element_text(size = 24),
  legend.text = element_text(size = 20),
  strip.text = element_text(size = 24),
  # panel.background = element_rect(fill = "white"),
  # panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
)

dark_theme <- main_theme +
  theme(
    plot.background  = element_rect(fill = "#060606", color = NA),
    panel.background = element_rect(fill = "#060606"),
    panel.border     = element_rect(color = "#282828", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "#282828"),
    panel.grid.minor = element_line(color = "#1a1a1a"),
    axis.ticks       = element_line(color = "#adafae"),
    axis.text        = element_text(color = "#adafae"),
    axis.title       = element_text(color = "#adafae"),
    plot.title       = element_text(color = "#ffffff"),
    strip.background = element_rect(fill = "#282828"),
    strip.text       = element_text(color = "#adafae"),
    legend.background = element_rect(fill = "#060606", color = NA),
    legend.key        = element_rect(fill = "#060606"),
    legend.title      = element_text(color = "#adafae"),
    legend.text       = element_text(color = "#adafae"),
  )

color_palette <- c(
  "#2a9fd6", "#e83e8c", "#93c", "#f80", "#77b300", "#20c997", "#c00", "#060606"
)
