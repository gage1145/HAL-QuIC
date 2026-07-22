library(tidyverse)
library(quicR)
library(janitor)
library(zoo)
library(modelr)
library(ggtext)
library(latex2exp)
library(ggpubr)

source("R/normalize.R")
source("R/fit_model.R")
source("R/estimate_params.R")
source("R/main_theme.R")

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

df_ <- df_raw %>%
  nest(.by = all_of(group_list)) %>%
  left_join(
    calculate_metrics(
      df_raw, group_list, threshold = 3,
      time_col = "time", ttt_values = "norm", auc_values = "norm",
      norm_col = "norm", deriv_col = "deriv"
    ),
    by = group_list
  ) %>%
  filter(MPR > 3) %>%
  estimate_params() %>%
  na.omit()


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
    rse = map_dbl(model, sigma),
    aic = map_dbl(model, AIC),
    bic = map_dbl(model, BIC),
    pseudo_r2 = map_dbl(data, ~ cor(.x$pred, .x$norm)^2),
  ) %>%
  unnest_wider(coefficients) %>%
  mutate(
    Kapp = 1 / b1,
    tlagprime = b1 - a1 / 2,
    intercept = (S1 / 2) - a1 * b1,
    # tlag = -intercept / a1,
    tlag = b1 - (S1 / (2 * a1))
  )

# Residuals
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

# Worst fits
df_worst <- df_mod %>%
  arrange(pseudo_r2) %>%
  slice_head(n = 12)

df_worst %>%
  unnest(data) %>%
  pivot_longer(c(pred, growth, decay), names_to = "series") %>%
  mutate(
    series = factor(
      series, 
      levels = c("growth", "decay", "pred"),
      labels = c("Growth", "Decay", "Prediction")
    ),
  ) %>%
  ggplot(aes(time, value, color = series)) +
  geom_point(aes(y = norm), size = 0.5, color = "black") +
  geom_line(linewidth = 1) +
  geom_vline(aes(xintercept = time_to_growth_max), linetype = "dashed", linewidth = 1) +
  facet_wrap(vars(rxn, wells), ncol = 4) +
  scale_x_continuous(breaks = seq(0, 70, by = 8)) +
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



# Highlighting tlag
df_sample <- df_mod %>%
  arrange(bic) %>%
  filter(b1 < 48) %>%
  slice_head(n = 12) %>%
  arrange(tlag) %>%
  mutate(
    # dilutions = as.factor(round(dilutions, 2)),
    label = paste0("t<sub>lag</sub> = ", round(tlag, 2), "h"),
    label_x = ifelse(tlag > 36, tlag - 2, tlag + 2),
    hjust = ifelse(tlag > 36, 1, 0),
    y = map2_dbl(model, tlag, ~ unlist(predict(.x, newdata = tibble(time = .y)))),
    label_y = ifelse(tlag > 36, y + 0.5, y),
    facet = 1:12
  ) 

df_sample %>% 
  # arrange(b1) %>%
  unnest(data) %>%
  mutate(
    tlag_curve = a1 * time + intercept,
    tlag_curve = ifelse(tlag_curve < min(decay) | tlag_curve > max(norm), NA, tlag_curve),
  ) %>%
  pivot_longer(c(pred, growth, decay), names_to = "series") %>%
  mutate(
    series = factor(series, levels = c("growth", "decay", "pred"), labels = c("Growth", "Decay", "Prediction")),
  ) %>%
  ggplot(aes(time, value, color = series)) +
  geom_point(aes(y = norm), size = 1, color = "black") +
  geom_line(linewidth = 1.5) +
  geom_line(aes(y = tlag_curve), linetype = "dashed", linewidth = 1.2, color = "black") +
  geom_segment(aes(y = 0, x = tlag, yend = y), data = df_sample, linetype = "dashed", linewidth = 1.2, inherit.aes = FALSE) +
  geom_segment(aes(y = y, x = 0, xend = label_x), data = df_sample, linetype = "dashed", linewidth = 1.2, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linewidth = 1) +
  geom_richtext(aes(label = label, x = label_x, y = label_y, hjust = hjust), data = df_sample, size = 8, label.size = NA, fill = NA, inherit.aes = FALSE) +
  facet_wrap(vars(facet), ncol = 4) +
  scale_x_continuous(breaks = seq(0, 70, by = 8)) +
  labs(x = "Time (hr)", y = "Normalized Fluorescence") +
  main_theme +
  theme(
    strip.text = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "inside",
    legend.position.inside = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.direction = "horizontal",
  )
ggsave("figures/kinetics/sample.png", width = 18, height = 12)
