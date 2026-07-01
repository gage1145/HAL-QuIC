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

group_list <- c("sample", "wells", "dilutions", "assay", "reaction", "mortem", "sample_type", "animal")

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

df_ <- raw_file %>%
  read_parquet() %>%
  mutate(across(all_of(group_list), as.factor)) %>%
  select(-c(norm, deriv)) %>%
  filter(time <= 72) %>%
  norm_n_der("time", "rfu", 8, group_list) %>%
  filter(norm > 4) %>%
  group_by(across(all_of(group_list))) %>%
  na.omit()

df_temp <- df_ %>%
  summarize(
    max_time             = max(time, na.rm=TRUE),
    max_deriv            = max(deriv, na.rm=TRUE),
    growth_scale         = max_deriv,
    peak_norm            = norm[which.min(deriv2)],
    time_to_growth_max   = time[which(norm == peak_norm)[1]],
    time_to_growth_mid   = time[which.max(deriv)],
    max_equillibrium     = max(norm[which(time >= time_to_growth_max)], na.rm=TRUE),
    min_equillibrium     = min(norm[which(time >= time_to_growth_max)], na.rm=TRUE),
    max_decay            = max_equillibrium - peak_norm,
    min_decay            = min_equillibrium - peak_norm,
    equillibrium         = ifelse(abs(min_decay) > max_decay, min_equillibrium, max_equillibrium),
    time_to_equillibrium = time[which(norm == equillibrium)][1],
    peak_decay           = equillibrium - peak_norm,
    time_to_decay        = time_to_equillibrium - time_to_growth_max,
    time_to_decay_mid    = time_to_growth_max + time_to_decay / 2,
    decay_slope          = peak_decay / time_to_decay,
    decay_slope          = replace_na(decay_slope, 0),
    decay_scale          = abs(decay_slope),
    .groups = "drop"
  )

df_ <- df_ %>%
  nest() %>%
  left_join(df_temp)

rm(df_temp)

fit_model <- function(data, 
                      peak_norm, 
                      time_to_growth_mid, 
                      growth_scale, 
                      peak_decay, 
                      time_to_decay_mid, 
                      decay_scale,
                      max_time,
                      ...) {
  form1 <- norm ~ (S1 / (1 + exp(a1 * (b1 - time))))
  form2 <- as.formula(paste(deparse(form1), "+ (S2 / (1 + exp(a2 * (b2 - time))))"))

  peak_scalar <- 2

  lower_S1 <- 0
  lower_a1 <- 0.01
  lower_b1 <- 0
  lower_S2 <- -peak_norm * peak_scalar
  lower_a2 <- 0
  lower_b2 <- time_to_growth_mid

  upper_S1 <- peak_norm * peak_scalar
  upper_a1 <- 20
  upper_b1 <- max_time
  upper_S2 <- peak_norm * peak_scalar
  upper_a2 <- 10
  upper_b2 <- max_time

  fit_single <- function() {
    single_mod <- NULL
    tryCatch(
      {
        single_mod <- nls(
          form1, data = data,
          start = list(
            S1 = peak_norm,
            a1 = peak_norm,
            b1 = time_to_growth_mid
          ),
          algorithm = "port",
          lower = c(S1 = lower_S1, a1 = lower_a1, b1 = lower_b1),
          upper = c(S1 = upper_S1, a1 = upper_a1, b1 = upper_b1)
        )
      },
      # silent = TRUE
    )
    return(coef(single_mod))
  }

  fit_double <- function(single_mod) {

    double_mod <- NULL

    try(
      {
        double_mod <- nls(
          form2, data = data,
          start = list(
            S1 = peak_norm,  a1 = growth_scale, b1 = time_to_growth_mid,
            S2 = peak_decay, a2 = decay_scale,  b2 = time_to_decay_mid
          ),
          algorithm = "port",
          lower = c(
            S1 = lower_S1, a1 = lower_a1, b1 = lower_b1, 
            S2 = lower_S2, a2 = lower_a2, b2 = lower_b2
          ),
          upper = c(
            S1 = upper_S1, a1 = upper_a1, b1 = upper_b1, 
            S2 = upper_S2, a2 = upper_a2, b2 = upper_b2
          )
        )
      },
      # silent = TRUE
    )

    if(is.null(double_mod) & !is.null(single_mod)) {
      try(
        {
          # coefs <- coef(single_mod)
          double_mod <- nls(
            form2, data = data,
            start = list(
              S1 = single_mod[1], a1 = single_mod[2], b1 = single_mod[3],
              S2 = peak_decay, a2 = decay_scale, b2 = time_to_decay_mid
            ),
            algorithm = "port",
            lower = c(
              S1 = lower_S1, a1 = lower_a1, b1 = lower_b1, 
              S2 = lower_S2, a2 = lower_a2, b2 = lower_b2
            ),
            upper = c(
              S1 = upper_S1, a1 = upper_a1, b1 = upper_b1, 
              S2 = upper_S2, a2 = upper_a2, b2 = upper_b2
            )
          )
        },
        # silent = TRUE
      )
    }
      
    if(is.null(double_mod)) return(single_mod)
    return(coef(double_mod))
  }

  mod <- NULL
  try(mod <- fit_single(), silent = TRUE)
  try(mod <- fit_double(mod), silent = TRUE)
  return(mod)
}

df_mod <- df_ %>%
  # ungroup() %>%
#   head(1) %>%
  slice_sample(n = 10) %>%
  mutate(
    model = pmap(., fit_model, .progress = TRUE)
  )

df_unmod <- df_mod %>%
  filter(map_lgl(model, is.null))

df_results <- df_mod %>%
  filter(map_lgl(model, ~ inherits(.x, "nls"))) %>%
  mutate(
    data = map2(model, data, ~ {
      cc <- coef(.x)
      .y %>%
        add_predictions(.x) %>%
        add_residuals(.x) %>%
        mutate(
          growth = cc[1] / (1 + exp((cc[3] - time) * cc[2])),
          decay  = cc[4] / (1 + exp((cc[6] - time) * cc[5]))
        )
    }),
    coefficients = map(model, coef)
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
  panel.border = element_rect(color = "black", fill = NA, size = 1)
)

df_unmod %>%
  ungroup() %>%
  unnest(data) %>%
  summarize(
    norm = mean(norm),
    .by = c(time, wells, assay, reaction)
  ) %>%
  ggplot(aes(time, norm)) +
  geom_line() +
  facet_wrap(vars(wells, assay, reaction)) +
  main_theme +
  theme(
    strip.text = element_blank(),
  )

# Samples with greatest deviation from model
df_long <- df_results %>%
  unnest(data)
  # ungroup() %>%
df_dev <- df_long %>%
group_by(across(all_of(group_list))) %>%
  summarize(
    dev = mean(resid, na.rm=TRUE),
    sd = sd(resid, na.rm=TRUE)
  )

overall_deviation <- mean(df_long$resid, na.rm=TRUE)
overall_sd <- sd(df_long$resid, na.rm=TRUE)
threshold <- abs(overall_deviation) + 0.15 * overall_sd

deviants <- df_dev %>%
  filter(abs(dev) > threshold) %>%
  arrange(desc(dev))

df_long %>%
  right_join(select(deviants, -c(dev, sd))) %>%
  pivot_longer(c(pred, resid, norm), names_to = "series") %>%
  mutate(series = factor(series, levels = c("norm", "pred", "resid"))) %>%
  ggplot(aes(time, value, color = series, linetype = series)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("black", "red", "blue")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
  facet_wrap(vars(reaction, assay, wells)) +
  labs(x = "Time", y = "Residuals", title = "Reactions with Highest Deviations from Model") +
  main_theme +
  theme(
    strip.text = element_blank(),
    legend.title = element_blank(),
  )
ggsave("figures/deviants.png", width = 16, height = 12)


# Residual Visualizations ------------------------------------------------


# Plot histogram of residuals

res_hist <- df_long %>%
  ggplot(aes(resid)) +
  geom_histogram(bins = 200, fill="black") +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Residuals", y = "Count", title = "Histogram of Residuals") +
  main_theme

# QQ Plot
qqplot <- df_long %>%
  ggplot(aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(-4, 4)) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "Normal Q-Q Plot") +
  main_theme

# Residuals over time
res_time <- df_long %>%
  na.omit() %>%
  mutate(diff = pred - norm) %>%
  summarize(
    .by = c(time, assay),
    mean = mean(diff, na.rm = TRUE),
    sd = sd(diff, na.rm = TRUE)
  ) %>%
  ggplot(aes(time, mean, color = assay, fill = assay)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2) +
  labs(
    x = "Time", y = "Residuals", title = "Residuals over Time",
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


# Example Fits -----------------------------------------------------------


df_results %>%
#   slice_sample(n=12) %>%
  arrange(peak_decay) %>%
  head(12) %>%
  mutate(
    across(c(S1,a1,b1,S2,a2,b2), ~ signif(., 2)),
    # label = TeX(sprintf(r"($f(t)=\frac{%s}{1+e^{%s(%s - t)}} + \frac{%s}{1+e^{%s(%s - t)}}$)", S1,a1,b1,S2,a2,b2), output = "character"),
  ) %>%
  unnest(data) %>%
  ggplot(aes(time)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(y=norm), size=0.1, color="black") +
  geom_line(aes(y=pred), linewidth=1.2, color="darkred", linetype="dashed") +
  geom_line(aes(y=growth), linewidth=1.2, color="darkgreen", linetype="dashed") +
  geom_line(aes(y=decay), linewidth=1.2, color="darkorange", linetype="dashed") +
  facet_wrap(vars(reaction, wells)) +
#   geom_text(aes(label = label), x = 0, y = -3, hjust = 0, inherit.aes = FALSE, size = 4, parse = TRUE) +
  labs(x = "Time (hr)", y = "Normalized Fluorescence", title = "Example Fits") +
  main_theme +
  theme(
    strip.text = element_blank(),
  )

ggsave("figures/example_fits.png", width = 20, height = 12)


# Save Results to Parquet ------------------------------------------------

write_parquet(select(df_results, -model), "data/results.parquet")
