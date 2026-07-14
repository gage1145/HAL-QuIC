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



raw_file <- "data/data.parquet"

group_list <- c("sample", "wells", "dilutions", "assay", "reaction", "mortem", "sample_type", "animal")

norm_n_der <- function(df, x, y, norm_point, groups, window=3, smooth=10, zero=TRUE) {
  df %>%
    mutate(
      norm   = rollmean(!!sym(y), smooth, na.pad=TRUE),
      norm   = norm / norm[norm_point] - ifelse(zero, 1, 0),
      deriv  = (lead(norm, window) - lag(norm, window)) / (lead(!!sym(x), window) - lag(!!sym(x), window)),
      deriv2 = (lead(deriv, window) - lag(deriv, window)) / (lead(!!sym(x), window) - lag(!!sym(x), window)),
      .by = all_of(groups)
    )
}

df_ <- raw_file %>%
  read_parquet() %>%
  mutate(across(all_of(group_list), as.factor)) %>%
  select(-c(norm, deriv)) %>%
  filter(time <= 72, sample == "P") %>%
  norm_n_der("time", "rfu", 8, group_list) %>%
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
  )
  # filter(peak_norm > 4)

df_combined <- df_ %>%
  nest(.by = group_list) %>%
  right_join(df_temp)

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

  fit_timeout   <- 10  # seconds per nls() call
  fit_control   <- nls.control(maxiter = 100, warnOnly = TRUE)
  nls_bounded <- function(...) {
    setTimeLimit(elapsed = fit_timeout, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
    nls(..., control = fit_control)
  }

  peak_scalar <- 2
  is_negative_decay <- peak_decay < 0

  lower_S1 <- peak_norm
  lower_a1 <- 0.1
  lower_b1 <- 0
  lower_S2 <- ifelse(is_negative_decay, -peak_norm * peak_scalar, 0)
  lower_a2 <- 0
  lower_b2 <- -time_to_decay_mid

  upper_S1 <- peak_norm * peak_scalar
  upper_a1 <- 20
  upper_b1 <- max_time
  upper_S2 <- ifelse(is_negative_decay, 0, peak_norm * peak_scalar)
  upper_a2 <- 10
  upper_b2 <- max_time

  fit_single <- function(...) {
    single_mod <- NULL
    single_mod <- nls_bounded(
      form1, data = data,
      start = list(
        S1 = peak_norm,
        a1 = peak_norm,
        b1 = time_to_growth_mid
      ),
      algorithm = "port",
      lower = c(S1 = lower_S1, a1 = lower_a1, b1 = lower_b1),
      upper = c(S1 = upper_S1, a1 = upper_a1, b1 = upper_b1)
    ) %>%
      try(silent = TRUE)

    return(single_mod)
  }

  fit_double <- function(...) {

    mod <- NULL
    mod <- nls_bounded(
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
    ) %>%
      try()

    if(is.null(mod)) {
      mod <- try(fit_single())
      
      if(is.null(mod)) return(NULL)
      
      coefs <- coef(mod)
      S1 <- coefs[1]
      a1 <- coefs[2]
      b1 <- coefs[3]

      mod <- nls_bounded(
        form2, data = data,
        start = list(
          S1 = S1, a1 = a1, b1 = b1,
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
      ) %>%
        try()
    }
      
    return(mod)
  }

  fit_double()
}

df_mod <- df_combined %>%
  mutate(model = pmap(., fit_model, .progress = TRUE))

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
          growth = cc[1] / (1 + exp(cc[2] * (cc[3] - time))),
          decay  = cc[4] / (1 + exp(cc[5] * (cc[6] - time)))
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
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
)

df_unmod %>%
  ungroup() %>%
  unnest(data) %>%
  summarize(
    norm = mean(norm),
    .by = c(time, wells, assay, reaction)
  ) %>%
  {
    ggplot(aes(time, norm)) +
    geom_line() +
    facet_wrap(vars(wells, assay, reaction)) +
    main_theme +
    theme(
      strip.text = element_blank(),
    ) 
  } %>%
  try(silent = TRUE)

# Samples with greatest deviation from model
df_long <- df_results %>%
  unnest(data)

overall_deviation <- mean(df_long$resid, na.rm=TRUE)
overall_sd <- sd(df_long$resid, na.rm=TRUE)
threshold_scalar <- 0.2
threshold <- abs(overall_deviation) + threshold_scalar * overall_sd

deviants <- df_long %>%
  summarize(
    dev = mean(resid, na.rm=TRUE),
    sd = sd(resid, na.rm=TRUE),
    .by = group_list
  ) %>%
  filter(abs(dev) > threshold) %>%
  arrange(desc(dev))

df_long %>%
  inner_join(select(deviants, -c(dev, sd))) %>%
  pivot_longer(c(pred, growth, decay, norm), names_to = "series") %>%
  mutate(series = factor(series, levels = c("norm", "pred", "growth", "decay"))) %>%
  ggplot(aes(time, value, color = series, linetype = series, linewidth = series)) +
  geom_line() +
  # geom_vline(aes(xintercept = time_to_growth_max), linetype = "dashed", linewidth=0.5) +
  # geom_vline(aes(xintercept = time_to_decay_mid), linetype = "dashed", color = "blue") +
  # geom_vline(aes(xintercept = time_to_equillibrium), linetype = "dashed", color = "orange") +
  scale_color_manual(values = c("black", "red", "orange", "blue")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed", "dashed")) +
  scale_linewidth_manual(values = c(1, 1.5, 1.5, 1.5)) +
  facet_wrap(vars(reaction, assay, wells)) +
  labs(x = "Time") +
  main_theme +
  theme(
    axis.title.y = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0, .95),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.direction = "horizontal",
  )
ggsave("figures/deviants.png", width = 16, height = 12)


# Residual Visualizations ------------------------------------------------


# Plot histogram of residuals

res_hist <- df_long %>%
  ggplot(aes(resid, fill = assay)) +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity") +
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

# QQ Plot
qqplot <- df_long %>%
  ggplot(aes(sample = resid, color = assay)) +
  geom_qq() +
  geom_qq_line() +
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(-4, 4)) +
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
#   filter(peak_norm > 4) %>%
  arrange(decay_slope) %>%
  head(24) %>%
  mutate(across(c(S1,a1,b1,S2,a2,b2), ~ signif(., 2))) %>%
  unnest(data) %>%
  ggplot(aes(time)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(aes(y=norm), linewidth=1, color="black") +
  geom_line(aes(y=pred), linewidth=1.5, color="red", linetype="dashed") +
  geom_line(aes(y=growth), linewidth=1.5, color="orange", linetype="dashed") +
  geom_line(aes(y=decay), linewidth=1.5, color="blue", linetype="dashed") +
  facet_wrap(vars(reaction, wells), ncol = 6, scales="free_y") +
  # geom_text(aes(label = label), x = 0, y = 30, hjust = 0, inherit.aes = FALSE, size = 4, parse = TRUE) +
  labs(x = "Time (hr)", y = "Normalized Fluorescence") +
  main_theme +
  theme(
    strip.text = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
  )

ggsave("figures/example_fits.png", width = 20, height = 12)


# Save Results to Parquet ------------------------------------------------

write_parquet(select(df_results, -model), "data/results.parquet")
