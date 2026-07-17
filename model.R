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

source("R/normalize.R")
source("R/fit_model.R")


single_only <- FALSE
weighted <- TRUE

# Noise model, sd^2 = s0^2 + (k * mu^theta)^2, estimated from the unweighted residuals. 
# Drives the IRLS weights in fit_model() and standardizes the residuals when ranking fits below.
irls_s0 <- 0.072
irls_k <- 0.0094
irls_theta <- 1

raw_file <- "data/data.parquet"

group_list <- c("sample", "wells", "dilutions", "assay", "reaction", "mortem", "sample_type", "animal")

df_ <- raw_file %>%
  read_parquet() %>%
  mutate(across(all_of(group_list), as.factor)) %>%
  select(-c(norm, deriv)) %>%
  filter(
    time <= 72,
    # sample == "P",
    assay == "RT-QuIC"
  ) %>%
  normalize("time", "rfu", 8, group_list) %>%
  na.omit()

df_temp <- df_ %>%
  summarize(
    max_time             = max(time, na.rm=TRUE),
    growth_scale         = max(deriv, na.rm=TRUE),
    peak_norm            = norm[which.min(deriv2)],
    time_to_growth_max   = time[norm == peak_norm][1],
    time_to_growth_mid   = time[which.max(deriv)],
    max_equillibrium     = max(norm[time > time_to_growth_max], na.rm=TRUE),
    min_equillibrium     = min(norm[time > time_to_growth_max], na.rm=TRUE),
    max_decay            = max_equillibrium - peak_norm,
    min_decay            = min_equillibrium - peak_norm,
    equillibrium         = ifelse(abs(min_decay) > max_decay, min_equillibrium, max_equillibrium),
    time_to_equillibrium = time[norm == equillibrium][1],
    peak_decay           = equillibrium - peak_norm,
    time_to_decay        = time_to_equillibrium - time_to_growth_max,
    time_to_decay_mid    = time_to_growth_max + time_to_decay / 2,
    decay_slope          = replace_na(peak_decay / time_to_decay, 0),
    decay_scale          = abs(decay_slope),
    .by = group_list
  ) %>%
  filter(
    peak_norm >= 3,
    time_to_growth_max < 50
  )

df_combined <- df_ %>%
  nest(.by = group_list) %>%
  right_join(df_temp)

rm(df_temp)

df_mod <- df_combined %>%
  mutate(model = pmap(
    ., fit_model,
    single_only = single_only, weighted = weighted,
    irls_s0 = irls_s0, irls_k = irls_k, irls_theta = irls_theta, .progress = TRUE
  ))

df_results <- df_mod %>%
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


# Figures ----------------------------------------------------------------


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

df_long <- df_results %>%
  unnest(data)

# Residual Visualizations ------------------------------------------------


# Residual Histogram
# res_hist <- df_long %>%
#   ggplot(aes(resid)) +
#   geom_histogram(bins = 100, position = "identity", fill = "blue") +
#   scale_x_continuous(limits = c(-1, 1)) +
#   labs(x = "Residuals", y = "Count", title = "Histogram of Residuals") +
#   main_theme +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "inside",
#     legend.position.inside = c(0.1, .95),
#     legend.justification = c(0, 1),
#     legend.background = element_blank(),
#     legend.direction = "horizontal",
#   )

# # QQ Plot
# qqplot <- df_long %>%
#   ggplot(aes(sample = resid)) +
#   geom_qq(color = "blue") +
#   geom_qq_line() +
#   scale_x_continuous(limits = c(-4, 4)) +
#   scale_y_continuous(limits = c(-4, 4)) +
#   labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "Normal Q-Q Plot") +
#   main_theme +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "inside",
#     legend.position.inside = c(0.1, .95),
#     legend.justification = c(0, 1),
#     legend.background = element_blank(),
#     legend.direction = "horizontal",
#   )

# # Residuals over time
# res_time <- df_long %>%
#   summarize(
#     .by = c(time, assay),
#     mean = mean(resid, na.rm = TRUE),
#     sd = sd(resid, na.rm = TRUE)
#   ) %>%
#   ggplot(aes(time, mean)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_line(color = "blue") +
#   geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), color = "blue", fill = "blue", alpha = 0.2) +
#   labs(
#     x = "Time", y = "Residuals", title = "Residuals over Time",
#   ) +
#   main_theme +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "inside",
#     legend.position.inside = c(0.1, .95),
#     legend.justification = c(0, 1),
#     legend.background = element_blank(),
#     legend.direction = "horizontal",
#   )

# (res_hist | qqplot) / res_time
# ggsave("figures/residual_vis.png", width = 16, height = 12)


# Worst Fits -------------------------------------------------------------


n_examples <- 12

fit_quality <- df_results %>%
  mutate(
    converged = map_lgl(model, ~ isTRUE(.x$convInfo$isConv)),
    srmse = map_dbl(data, ~ {
      sd_hat <- sqrt(irls_s0^2 + (irls_k * abs(.x$pred)^irls_theta)^2)
      sqrt(mean((.x$resid / sd_hat)^2, na.rm = TRUE))
    }),
    # Worst single point: catches localized misfit.
    max_sres = map_dbl(data, ~ {
      sd_hat <- sqrt(irls_s0^2 + (irls_k * abs(.x$pred)^irls_theta)^2)
      max(abs(.x$resid / sd_hat), na.rm = TRUE)
    })
  )

deviants <- fit_quality %>%
  arrange(desc(srmse)) %>%
  head(n_examples) %>%
  select(all_of(group_list), srmse, max_sres, converged) %>%
  left_join(df_long, by = group_list)

df_pred <- deviants %>%
  select(group_list, time, pred, growth, decay) %>%
  pivot_longer(c(pred, growth, decay), names_to = "series") %>%
  mutate(series = factor(series, levels = c("pred", "growth", "decay")))

df_ci <- deviants %>%
  select(group_list, time, norm, lower, upper)

deviants %>%
  ggplot(aes(time, norm)) +
  geom_point(size = 0.5) +
  geom_line(aes(y = value, color = series), data = df_pred, linewidth = 1.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), data = df_ci, alpha = 0.2) +
  facet_wrap(vars(reaction, assay, wells)) +
  labs(x = "Time") +
  main_theme +
  theme(
    axis.title.y = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.position.inside = c(0, .95),
    legend.background = element_blank(),
    legend.direction = "horizontal",
  )
ggsave("figures/deviants.png", width = 8, height = 6)


# Best Fits -----------------------------------------------------------


fit_quality %>%
  arrange(srmse) %>%
  head(n_examples) %>%
  select(all_of(group_list), srmse, max_sres, converged) %>%
  left_join(df_long, by = group_list) %>%
  mutate(across(matches("$[A-z]{1}\\d^"), ~ signif(., 2))) %>%
  pivot_longer(c(pred, growth, decay), names_to = "series") %>%
  mutate(series = factor(series, levels = c("pred", "growth", "decay"))) %>%
  ggplot(aes(time)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(y = norm), size = 0.5, color = "black") +
  geom_line(aes(y = value, color = series), linewidth = 1.5, linetype = "dashed") +
  facet_wrap(vars(reaction, wells)) +
  labs(x = "Time (hr)", y = "Normalized Fluorescence") +
  main_theme +
  theme(
    axis.title.y = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.position.inside = c(0, .95),
    legend.background = element_blank(),
    legend.direction = "horizontal",
  )

ggsave("figures/best_fits.png", width = 8, height = 6)


# Save Results to Parquet ------------------------------------------------


write_parquet(select(fit_quality, -model), "data/results.parquet")
