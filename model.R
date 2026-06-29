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
  select(-c(norm, deriv)) %>%
  filter(time != 0)

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

fit_model <- function(data, 
                      peak_norm, 
                      time_to_growth_mid, 
                      growth_scale, 
                      peak_decay, 
                      time_to_decay_mid, 
                      decay_scale,
                      max_time,
                      ...) {
  form1 <- norm ~ (a / (1 + exp(c * (b - time))))
  form2 <- as.formula(paste(deparse(form1), "+ (d / (1 + exp(f * (e - time))))"))

  peak_scalar <- 2
  time_scalar <- 0.9

  lower_a <- 0
  lower_b <- 0
  lower_c <- 0.01
  lower_d <- -peak_norm * peak_scalar
  lower_e <- time_to_growth_mid
  lower_f <- 0

  upper_a <- peak_norm * peak_scalar
  upper_b <- max_time * time_scalar
  upper_c <- 10
  upper_d <- peak_norm * peak_scalar
  upper_e <- max_time * time_scalar
  upper_f <- 5

  fit_single <- function() {
    single_mod <- NULL
    tryCatch(
      {
        single_mod <- nls(
          form1, data = data,
          start = list(
            a = peak_norm,
            b = time_to_growth_mid,
            c = growth_scale
          ),
          algorithm = "port",
          lower = c(a = lower_a, b = lower_b, c = lower_c),
          upper = c(a = upper_a, b = upper_b, c = upper_c)
        )
      },
      error = function(e) {
        sprintf("Unable to fit single modlel: %s", e)
      }
    )
    return(single_mod)
  }

  fit_double <- function(single_mod) {

    double_mod <- NULL

    try({
      double_mod <- nls(
        form2, data = data,
        start = list(
          a = peak_norm,  b = time_to_growth_mid, c = growth_scale,
          d = peak_decay, e = time_to_decay_mid,  f = decay_scale
        ),
        algorithm = "port",
        lower = c(
          a = lower_a, b = lower_b, c = lower_c, 
          d = lower_d, e = lower_e, f = lower_f
        ),
        upper = c(
          a = upper_a, b = upper_b, c = upper_c, 
          d = upper_d, e = upper_e, f = upper_f
        )
      )
    })

    if(is.null(double_mod) & !is.null(single_mod)) {
      try({
        coefs <- coef(single_mod)
        double_mod <- nls(
          form2, data = data,
          start = list(
            a = coefs[1],
            b = coefs[2],
            c = coefs[3],
            d = peak_decay,
            e = time_to_decay_mid,
            f = decay_scale
          ),
          algorithm = "port",
          lower = c(
            a = lower_a, b = lower_b, c = lower_c, 
            d = lower_d, e = lower_e, f = lower_f
          ),
          upper = c(
            a = upper_a, b = upper_b, c = upper_c, 
            d = upper_d, e = upper_e, f = upper_f
          )
        )
      })
    }
      
    if(is.null(double_mod)) return(single_mod)
    return(double_mod)
  }

  mod <- NULL
  try(mod <- fit_single(), silent = TRUE)
  try(mod <- fit_double(mod), silent = TRUE)
  return(mod)
}

df_mod <- df_test %>%
  group_by(across(all_of(group_list))) %>%
  nest() %>%
  ungroup() %>%
  left_join(df_test_sum) %>%
  mutate(
    model = pmap(., fit_model)
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
          growth = cc[1] / (1 + exp((cc[2] - time) * cc[3])),
          decay  = cc[4] / (1 + exp((cc[5] - time) * cc[6]))
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

# df_unmod %>%
#   ungroup() %>%
#   unnest(data) %>%
#   summarize(
#     norm = mean(norm),
#     .by = c(time, wells, assay, reaction)
#   ) %>%
#   ggplot(aes(time, norm)) +
#   geom_line() +
#   facet_wrap(vars(wells, assay, reaction)) +
#   main_theme +
#   theme(
#     strip.text = element_blank(),
#   )

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
  main_theme
res_time

(res_hist | qqplot) / res_time
ggsave("figures/residual_vis.png", width = 16, height = 12)


# Example Fits -----------------------------------------------------------


df_results %>%
#   slice_sample(n=12) %>%
  arrange(decay_slope) %>%
  head(12) %>%
  mutate(
    across(c(a,b,c,d,e,f), ~ signif(., 2)),
    label = TeX(sprintf(r"($f(t)=\frac{%s}{1+e^{%s(%s - t)}} + \frac{%s}{1+e^{%s(%s - t)}}$)", a,b,c,d,e,f), output = "character"),
  ) %>%
  unnest(data) %>%
  ggplot(aes(time)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(aes(y=norm), linewidth=0.5, color="blue") +
  geom_line(aes(y=pred), linewidth=1.2, color="darkred", linetype="dashed") +
  geom_line(aes(y=growth), linewidth=1.2, color="darkgreen", linetype="dashed") +
  geom_line(aes(y=decay), linewidth=1.2, color="darkorange", linetype="dashed") +
  facet_wrap(vars(reaction, wells)) +
  geom_text(aes(label = label), x = 30, y = 2, hjust = 0, inherit.aes = FALSE, size = 4, parse = TRUE) +
  labs(x = "Time (hr)", y = "Normalized Fluorescence", title = "Example Fits") +
  main_theme +
  theme(
    strip.text = element_blank(),
  )

ggsave("figures/example_fits.png", width = 20, height = 12)
