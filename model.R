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


single_only <- FALSE
weighted <- TRUE

# Noise model, sd^2 = s0^2 + (k * mu^theta)^2, estimated from the unweighted
# residuals. Drives the IRLS weights in fit_model() and standardizes the
# residuals when ranking fits below -- keep them a single source of truth.
irls_s0 <- 0.072
irls_k <- 0.0094
irls_theta <- 1

raw_file <- "data/data.parquet"

group_list <- c("sample", "wells", "dilutions", "assay", "reaction", "mortem", "sample_type", "animal")

norm_n_der <- function(df, x, y, norm_point, groups, window = 3, smooth = 10, zero = TRUE) {
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
  filter(
    time <= 72,
    sample == "P",
    assay == "RT-QuIC"
  ) %>%
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
                      ...,
                      single_only = FALSE,
                      algorithm = "port",
                      weighted = FALSE,
                      irls_iter = 3,
                      irls_s0 = 0.072,
                      irls_k = 0.0094,
                      irls_theta = 1) {
  form1 <- norm ~ (S1 / (1 + exp(a1 * (b1 - time))))
  form2 <- as.formula(paste(deparse(form1), "+ (S2 / (1 + exp(a2 * (b2 - time))))"))

  data$w_vector <- 1

  fit_timeout <- 10 # seconds per nls() call
  fit_control <- nls.control(maxiter = 1000, warnOnly = TRUE)
  nls_bounded <- function(...) {
    setTimeLimit(elapsed = fit_timeout, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
    nls(..., control = fit_control, algorithm = algorithm, weights = w_vector)
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

  fit_single <- function(start = NULL, ...) {
    if (is.null(start)) {
      start <- list(
        S1 = peak_norm,
        a1 = peak_norm,
        b1 = time_to_growth_mid
      )
    }

    single_mod <- NULL
    single_mod <- nls_bounded(
      form1,
      data = data,
      start = start,
      lower = c(S1 = lower_S1, a1 = lower_a1, b1 = lower_b1),
      upper = c(S1 = upper_S1, a1 = upper_a1, b1 = upper_b1)
    ) %>%
      try(silent = TRUE)

    return(single_mod)
  }

  # Iteratively reweighted least squares against the variance model
  irls <- function(fit_fun) {
    mod <- fit_fun()
    if (!weighted || !inherits(mod, "nls")) {
      return(mod)
    }

    for (i in seq_len(irls_iter)) {
      w <- 1 / (irls_s0^2 + (irls_k * abs(predict(mod))^irls_theta)^2)
      if (!all(is.finite(w)) || any(w <= 0)) break
      data$w_vector <<- w

      new_mod <- fit_fun(start = as.list(coef(mod)))
      if (!inherits(new_mod, "nls")) break # refit failed; keep last good fit
      mod <- new_mod
    }

    return(mod)
  }

  if (single_only) {
    return(irls(fit_single))
  }

  fit_double <- function(start = NULL, ...) {
    if (is.null(start)) {
      start <- list(
        S1 = peak_norm,  a1 = growth_scale, b1 = time_to_growth_mid,
        S2 = peak_decay, a2 = decay_scale,  b2 = time_to_decay_mid
      )
    }

    mod <- NULL
    mod <- nls_bounded(
      form2,
      data = data,
      start = start,
      lower = c(
        S1 = lower_S1, a1 = lower_a1, b1 = lower_b1,
        S2 = lower_S2, a2 = lower_a2, b2 = lower_b2
      ),
      upper = c(
        S1 = upper_S1, a1 = upper_a1, b1 = upper_b1,
        S2 = upper_S2, a2 = upper_a2, b2 = upper_b2
      )
    ) %>%
      try(silent = TRUE)

    if (is.null(mod)) {
      mod <- try(fit_single())

      if (is.null(mod)) {
        return(NULL)
      }

      coefs <- coef(mod)
      S1 <- coefs[1]
      a1 <- coefs[2]
      b1 <- coefs[3]

      mod <- nls_bounded(
        form2,
        data = data,
        start = list(
          S1 = S1, a1 = a1, b1 = b1,
          S2 = peak_decay, a2 = decay_scale, b2 = time_to_decay_mid
        ),
        lower = c(
          S1 = lower_S1, a1 = lower_a1, b1 = lower_b1,
          S2 = lower_S2, a2 = lower_a2, b2 = lower_b2
        ),
        upper = c(
          S1 = upper_S1, a1 = upper_a1, b1 = upper_b1,
          S2 = upper_S2, a2 = upper_a2, b2 = upper_b2
        )
      ) %>%
        try(silent = TRUE)
    }

    return(mod)
  }

  irls(fit_double)
}

df_mod <- df_combined %>%
  mutate(model = pmap(
    ., fit_model,
    single_only = single_only, weighted = weighted,
    irls_s0 = irls_s0, irls_k = irls_k, irls_theta = irls_theta, .progress = TRUE
  ))

df_unmod <- df_mod %>%
  filter(map_lgl(model, is.null))

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

df_long <- df_results %>%
  unnest(data)

# Residual Visualizations ------------------------------------------------


# Residual Histogram
res_hist <- df_long %>%
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

# QQ Plot
qqplot <- df_long %>%
  ggplot(aes(sample = resid)) +
  geom_qq(color = "blue") +
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
  summarize(
    .by = c(time, assay),
    mean = mean(resid, na.rm = TRUE),
    sd = sd(resid, na.rm = TRUE)
  ) %>%
  ggplot(aes(time, mean)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), color = "blue", fill = "blue", alpha = 0.2) +
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
res_time

(res_hist | qqplot) / res_time
ggsave("figures/residual_vis.png", width = 16, height = 12)


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
