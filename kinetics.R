library(tidyverse)
library(quicR)
library(janitor)
library(zoo)
library(modelr)
library(ggtext)

source("R/normalize.R")
source("R/fit_model.R")
source("R/estimate_params.R")

main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 20),
  legend.title = element_text(size = 24),
  legend.text = element_text(size = 20),
  strip.text = element_text(size = 24),
  panel.background = element_rect(fill = "white"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
)

files <- list.files("raw/kinetics", pattern = "*.xlsx", full.names = TRUE)
group_list <- c("sample", "wells", "dilutions", "rxn")

get_raw <- function(file) {
  file %>%
    get_quic() %>%
    mutate(
      rxn = str_extract(file, "2026\\d{4}_r\\d+_[A-Z]+_\\D+\\."),
      rxn = str_remove(rxn, "\\.")
    )
}

df_raw <- files %>%
  map_dfr(get_raw) %>%
  clean_names(replace = c("Sample IDs" = "sample")) %>%
  filter(time > 0) %>%
  normalize("time", "rfu", 8, group_list, 3) %>%
  mutate(
    dilutions = -log10(as.integer(dilutions)),
  ) %>%
  na.omit()

df_ <- 
  nest(.by = all_of(group_list)) %>%
  left_join(
    calculate_metrics(
      df_, group_list, threshold = 3,
      time_col = "time", ttt_values = "norm", auc_values = "norm",
      norm_col = "norm", deriv_col = "deriv"
    ),
    by = group_list
  ) %>%
  filter(MPR > 3) %>%
  estimate_params()

df_ %>%
  unnest(data) %>%
  # mutate(dilutions = as.factor(round(dilutions, 2))) %>%
  # summarize(norm = median(norm), .by = all_of(c("time", "sample", "dilutions", "rxn"))) %>%
  ggplot(aes(time, norm, color = dilutions, group = wells)) +
  geom_line(alpha = 0.5, linewidth = 1.5) +
  facet_grid(vars(rxn), scales = "free_y") +
  scale_color_gradient(low = "#0072B2", high = "#D55E00") +
  main_theme +
  theme(
    legend.position = "bottom"
  )



df_metrics %>%
  # summarize(RAF = mean(RAF), .by = c(sample, dilutions, rxn)) %>%
  ggplot(aes(dilutions, RAF)) +
  geom_point(size = 1.2) +
  stat_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.5)


# Modeling ---------------------------------------------------------------


df_mod <- df_ %>%
  mutate(model = pmap(
    ., fit_model,
    single_only = FALSE, weighted = TRUE, .progress = TRUE
  )) %>%
  filter(map_lgl(model, ~ inherits(.x, "nls"))) %>%
  mutate(
    data = map2(model, data, ~ {
      cc <- coef(.x)
      ci <- tryCatch(
        predFit(.x, newdata = .y, interval = "confidence", level = 0.95),
        error = function(e) NULL
      )
      .y %>%
        add_predictions(.x) %>%
        add_residuals(.x) %>%
        mutate(
          lower = if (is.null(ci)) NA_real_ else ci[, "lwr"],
          upper = if (is.null(ci)) NA_real_ else ci[, "upr"],
          growth = cc[1] / (1 + exp(cc[2] * (cc[3] - time))),
          decay  = cc[4] / (1 + exp(cc[5] * (cc[6] - time))),
        )
    }),
    coefficients = map(model, coef),
  ) %>%
  unnest_wider(coefficients)


df_mod %>%
  unnest(data) %>%
  ggplot(aes(resid)) +
  geom_histogram(bins = 100, position = "identity", fill = "blue") +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Residuals", y = "Count", title = "Histogram of Residuals") +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.1, .95),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.direction = "horizontal",
  )

df_sample <- df_mod %>%
  slice_sample(n = 12) %>%
  mutate(
    label = paste0("t<sub>lag</sub> = ", round(b1, 2), "h"),
    label_x = ifelse(b1 > 36, b1 - 2, b1 + 2),
    hjust = ifelse(b1 > 36, 1, 0),
    y = map2_dbl(model, b1, ~ unlist(predict(.x, newdata = tibble(time = .y)))),
    label_y = ifelse(b1 > 36, y + (peak_norm - y) / 2, y),
    facet = as.factor(1:12)
  ) 

df_sample %>% 
  arrange(b1) %>%
  unnest(data) %>%
  pivot_longer(c(pred, growth, decay), names_to = "series") %>%
  mutate(
    series = factor(series, levels = c("growth", "decay", "pred"), labels = c("Growth", "Decay", "Prediction")),
  ) %>%
  ggplot(aes(time, value, color = series)) +
  geom_point(aes(y = norm), size = 2, color = "black") +
  geom_line(linewidth = 1.5) +
  geom_vline(aes(xintercept = b1), data = df_sample, linetype = "dashed") +
  geom_segment(aes(y = y, x = 0, xend = label_x), data = df_sample, linetype = "dashed", inherit.aes = FALSE) +
  geom_richtext(aes(label = label, x = label_x, y = label_y, hjust = hjust), data = df_sample, size = 8, label.size = NA, fill = NA, inherit.aes = FALSE) +
  facet_wrap(vars(fct_inorder(facet)), ncol = 4) +
  labs(x = "Time (hr)", y = "Normalized Fluorescence") +
  main_theme +
  theme(
    strip.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.direction = "horizontal",
  )
ggsave("figures/kinetics/sample.png", width = 18, height = 12)
