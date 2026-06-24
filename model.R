library(tidyverse)
library(broom)
library(arrow)
library(magrittr)
library(modelr)
library(quicR)
library(zoo)
library(skimr)
library(ggpubr)
library(patchwork)
library(latex2exp)



raw_file <- "data/data.parquet"

group_list <- c("sample", "wells", "dilutions", "assay", "reaction")

df_ <- raw_file %>%
  read_parquet() %>%
  mutate(
    across(c(reaction, wells, dilutions, assay, mortem, sample_type, animal, sample), as.factor)
  ) %>%
  select(-c(norm, deriv))

df_p <- df_ %>%
  filter(sample %in% c("P"), time <= 72) 

norm_n_der <- function(df, x, y, norm_point, groups, window=3, smooth=10, zero=TRUE) {
  df %>%
    group_by(across(all_of(groups))) %>%
    mutate(
      norm   = rollmean(!!sym(y), smooth, na.pad=TRUE),
      norm   = norm / norm[norm_point] - ifelse(zero, 1, 0),
      deriv  = (lead(norm, window) - lag(norm, window)) / (lead(!!sym(x), window) - lag(!!sym(x), window)),
      deriv2 = (lead(deriv, window) - lag(deriv, window)) / (lead(!!sym(x), window) - lag(!!sym(x), window))
    )
}

df_test <- df_p %>%
  norm_n_der("time", "rfu", 8, group_list)

df_test_sum <- df_test %>%
  group_by(across(all_of(group_list))) %>%
  na.omit() %>%
  summarize(
    max_time             = max(time, na.rm=TRUE),
    max_deriv            = max(deriv, na.rm=TRUE),
    growth_scale         = max_deriv,
    peak_norm            = norm[which.min(deriv2)],
    time_to_growth_max   = time[which(norm == peak_norm)[1]],
    time_to_growth_mid   = time[which.max(deriv)],
    max_equillibrium     = max(norm[which(time > time_to_growth_max)], na.rm=TRUE),
    min_equillibrium     = min(norm[which(time > time_to_growth_max)], na.rm=TRUE),
    max_decay            = max_equillibrium - peak_norm,
    min_decay            = min_equillibrium - peak_norm,
    equillibrium         = ifelse(abs(min_decay) > max_decay, min_equillibrium, max_equillibrium),
    time_to_equillibrium = time[which(norm == equillibrium)][1],
    peak_decay           = equillibrium - peak_norm,
    time_to_decay        = time_to_equillibrium - time_to_growth_max,
    time_to_decay_mid    = time_to_growth_max + time_to_decay / 2,
    decay_slope          = peak_decay / time_to_decay,
    decay_scale          = abs(decay_slope),
    .groups = "drop"
  )



# form3 <- as.formula(paste(deparse(form2), "+ (g / (1 + exp(i * (h - time))))"))

fit_model <- function(df) {
  form1 <- norm ~ (a / (1 + exp(c * (b - time))))
  form2 <- as.formula(paste(deparse(form1), "+ (d / (1 + exp(f * (e - time))))"))

  fit_single <- function(df) {
    nls(
      form1, data = df,
      start = list(
        a = df$peak_norm[1],
        b = df$time_to_growth_mid[1],
        c = df$growth_scale[1]
      ),
      algorithm = "port",
      lower = c(a = 0, b = 0, c = 0),
      upper = c(a = df$peak_norm[1] * 1.5, b = df$max_time[1], c = 20)
    )
  }

  fit_double <- function(df) {
    nls(
      form2, data = df,
      start = list(
        a = df$peak_norm[1],
        b = df$time_to_growth_mid[1],
        c = df$growth_scale[1],
        d = df$peak_decay[1],
        e = df$time_to_decay_mid[1],
        f = df$decay_scale[1]
      ),
      algorithm = "port",
      lower = c(a = 0, b = 0, c = 0, d = -df$peak_norm[1] * 1.5, e = df$time_to_growth_mid[1], f = 0),
      upper = c(a = df$peak_norm[1] * 1.5, b = df$max_time[1], c = 20, d = df$peak_norm[1] * 1.5, e = df$max_time[1], f = 5)
    )
  }

  try(return(fit_double(df)), silent = TRUE)
  try(return(fit_single(df)), silent = FALSE)
  return(NULL)
}

df_mod <- df_test %>%
  left_join(df_test_sum) %>%
  group_by(across(all_of(group_list))) %>%
  nest() %>%
  mutate(
    model = map(data, fit_model)
  )

df_results <- df_mod %>%
  filter(map_lgl(model, ~ inherits(.x, "nls"))) %>%
  mutate(
    augmented = map2(model, data, ~ {
      cc <- coef(.x)
      .y %>%
        add_predictions(.x) %>%
        add_residuals(.x) %>%
        mutate(
          growth = cc[1] / (1 + exp((cc[2] - time) / cc[3])),
          decay  = cc[4] / (1 + exp((cc[5] - time) / cc[6]))
        )
    })
  ) %>%
  unnest(augmented)

df_sum <- df_mod %>%
  group_by(across(all_of(group_list))) %>%
  summarize(
    cc = map(model, coef)
  ) %>%
  unnest(cc) %>%
  mutate(coefficient = names(cc)) %>%
  pivot_wider(names_from = coefficient, values_from = cc) 

df_unmod <- df_mod %>%
  filter(map_lgl(model, is.null))

df_unmod %>%
  ungroup() %>%
  unnest(data) %>%
  summarize(
    norm = mean(norm),
    .by = c(time, wells, assay)
  ) %>%
  ggplot(aes(time, norm)) +
  geom_line() +
  facet_grid(vars(wells), vars(assay))


# Figures ----------------------------------------------------------------


# Samples with greatest deviation from model
df_dev <- df_results %>%
  # ungroup() %>%
  group_by(across(all_of(group_list))) %>%
  summarize(
    dev = mean(resid, na.rm=TRUE),
    sd = sd(resid, na.rm=TRUE)
  )

overall_deviation <- mean(df_results$resid, na.rm=TRUE)
overall_sd <- sd(df_results$resid, na.rm=TRUE)
threshold <- overall_deviation + overall_sd

deviants <- df_dev %>%
  filter(dev > threshold)

df_results %>%
  right_join(select(deviants, -c(dev, sd))) %>%
  pivot_longer(c(pred, resid, norm), names_to = "series") %>%
  mutate(series = factor(series, levels = c("norm", "pred", "resid"))) %>%
  ggplot(aes(time, value, color = series, linetype = series)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("black", "red", "blue")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
  facet_wrap(vars(reaction, assay)) +
  labs(x = "Time", y = "Residuals", title = "Reactions with Significant Deviations from Model") +
  main_theme +
  theme(
    strip.text = element_blank(),
    legend.title = element_blank(),
  )


# Residual Visualizations ------------------------------------------------


main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 20),
  legend.title = element_text(size = 24),
  legend.text = element_text(size = 20),
  strip.text = element_text(size = 24),
  panel.background = element_rect(fill = "white"),
  panel.border = element_rect(color = "black", fill = NA, size = 1)
)

# Plot histogram of residuals
res_hist <- df_results %>%
  ggplot(aes(resid)) +
  geom_histogram(bins = 400, fill="black") +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Residuals", y = "Count", title = "Histogram of Residuals") +
  main_theme

# QQ Plot
qqplot <- df_results %>%
  ggplot(aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "Normal Q-Q Plot") +
  main_theme

# Residuals over time
res_time <- df_results %>%
  ggplot(aes(time, resid)) +
  geom_bin2d(bins = 50, color = NULL) +
  scale_fill_gradient(low = "white", high="darkred") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, 4), expand=expansion()) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2)) +
  labs(x = "Time", y = "Residual", title = "Residuals over Time", fill="Count") +
  main_theme +
  theme(
    strip.text = element_blank(),
    legend.key.width = unit(6, "cm"),
    # panel.background = element_rect(fill = "black"),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

(res_hist | qqplot) / res_time
ggsave("figures/residual_vis.png", width = 16, height = 12)



# Example Fits -----------------------------------------------------------


df_sum <- df_sum %>%
  mutate(
    across(c(a,b,c,d,e,f), ~ signif(., 2)),
    label = TeX(sprintf("$f(t)=\\frac{%s}{(1+e^{%s(%s - t)}} + \\frac{%s}{(1+e^{%s(%s - t)}}$", a,b,c,d,e,f), output = "character"),
  )

df_sum <- df_sum %>%
  filter(reaction %in% unique(df_results$reaction)[1:10]) 

df_results %>%
  filter(reaction %in% unique(df_results$reaction)[1:10]) %>%
  ggplot(aes(time)) +
  geom_line(aes(y=norm), linewidth=0.5, color="blue") +
  geom_line(aes(y=pred), linewidth=1.2, color="darkred", linetype="dashed") +
  facet_grid(cols=vars(reaction), rows=vars(wells)) +
  geom_text(aes(label = label), data = df_sum, x = 0, y = 10, hjust = 0, inherit.aes = FALSE, size = 4, parse = TRUE) +
  scale_y_continuous(limits = c(0, 11)) +
  main_theme +
  theme(
    strip.text = element_blank(),
  )

ggsave("figures/example_fits.png", width = 28, height = 12)
