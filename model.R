library(tidyverse)
library(broom)
library(arrow)
library(magrittr)
library(modelr)
library(zoo)
library(skimr)
library(ggpubr)
library(patchwork)
library(latex2exp)
library(investr)
library(quicR)

source("R/normalize.R")
source("R/fit_model.R")
source("R/estimate_params.R")
source("R/main_theme.R")


single_only <- FALSE

raw_file <- "data/data.parquet"

group_list <- c("sample", "wells", "dilutions", "assay", "reaction", "mortem", "sample_type", "animal")

df_raw <- raw_file %>%
  read_parquet() %>%
  mutate(across(all_of(group_list), as.factor)) %>%
  filter(
    time <= 72,
    sample == "P",
    dilutions == -3,
    assay == "RT-QuIC"
  ) %>%
  normalize("time", "rfu", 8, group_list) %>%
  na.omit()

df_ <- df_raw %>%
  nest(.by = group_list) %>%
  left_join(
    calculate_metrics(
      df_raw, group_list, threshold = 3,
      time_col = "time", ttt_values = "norm", auc_values = "norm",
      norm_col = "norm", deriv_col = "deriv"
    ),
    by = group_list
  ) %>%
  estimate_params() %>%
  filter(MPR > 3 & time_to_growth_mid != max_time) %>%
  na.omit()

df_mod <- df_ %>%
  mutate(model = pmap(., fit_model, single_only = single_only, .progress = TRUE)) %>%
  filter(map_lgl(model, ~ inherits(.x, "nls"))) %>%
  mutate(
    data = map2(model, data, ~ {
      cc <- coef(.x)
      .y %>%
        add_predictions(.x) %>%
        add_residuals(.x) %>%
        mutate(
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
    tlag = b1 - a1 / 2
  )


# Residual Visualizations ------------------------------------------------


# Residual Histogram
res_hist <- df_mod %>%
  unnest(data) %>%
  # summarize(resid = mean(resid), .by = c(reaction, wells)) %>%
  ggplot(aes(resid)) +
  # geom_histogram(bins = 100, position = "identity", fill = "blue") +
  geom_density(fill = "blue") +
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
# res_hist

# QQ Plot
qqplot <- df_mod %>%
  unnest(data) %>%
  summarize(resid = mean(resid), .by = c(reaction, wells)) %>%
  ggplot(aes(sample = resid)) +
  geom_qq(color = "blue") +
  geom_qq_line() +
  scale_x_continuous(breaks = seq(-5, 5)) +
  # scale_y_continuous(breaks = seq(-5, 5)) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "Normal Q-Q Plot") +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.1, .95),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.direction = "horizontal",
  )
# qqplot

# Residuals over time
res_time <- df_mod %>%
  unnest(data) %>%
  summarize(
    .by = c(time),
    mean = mean(resid, na.rm = TRUE),
    sd = sd(resid, na.rm = TRUE)
  ) %>%
  ggplot(aes(time, mean)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), color = "blue", fill = "blue", alpha = 0.2) +
  labs(
    x = "Time", y = "Mean Residual", title = "Residuals over Time",
  ) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.1, .95),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.direction = "horizontal",
  )

(res_hist | qqplot) / res_time
ggsave("figures/residual_vis.png", width = 16, height = 12)


# Worst Fits -------------------------------------------------------------


n_examples <- 24

deviants <- df_mod %>%
  arrange(pseudo_r2) %>%
  slice_head(n = n_examples)

deviants %>%
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
  geom_vline(aes(xintercept = tlagprime), linetype = "dashed", linewidth = 1) +
  facet_wrap(vars(reaction, wells), ncol = 4) +
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
ggsave("figures/deviants.png", width = 18, height = 24)


# Best Fits -----------------------------------------------------------


df_mod %>%
  arrange(desc(pseudo_r2)) %>%
  filter(time_to_growth_max < 36) %>%
  head(n_examples) %>%
  mutate(across(matches("$[A-z]{1}\\d^"), ~ signif(., 2))) %>%
  unnest(data) %>%
  pivot_longer(c(pred, growth, decay), names_to = "series") %>%
  mutate(series = factor(series, levels = c("pred", "growth", "decay"))) %>%
  ggplot(aes(time)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(y = norm), size = 0.5, color = "black") +
  geom_line(aes(y = value, color = series), linewidth = 1) +
  facet_wrap(vars(reaction, wells), ncol = 4) +
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

ggsave("figures/best_fits.png", width = 18, height = 24)


# Save Results to Parquet ------------------------------------------------


write_parquet(select(df_mod, -model), "data/results.parquet")
