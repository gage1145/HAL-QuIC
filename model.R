library(tidyverse)
library(broom)
library(arrow)
library(magrittr)
library(modelr)
library(quicR)
library(zoo)
library(skimr)



raw_file <- "data/data.parquet"

group_list <- c("sample", "wells", "dilutions", "assay", "reaction")

df_ <- raw_file %>%
  read_parquet() %>%
  mutate(
    across(c(reaction, wells, dilutions, assay, mortem, sample_type, animal, sample), as.factor)
  ) %>%
  select(-c(norm, deriv))

df_p <- df_ %>%
  filter(sample %in% c("N", "P"), time <= 72) 

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
    growth_scale         = max_deriv / 4,
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
    decay_scale          = abs(decay_slope) * 100,
    .groups = "drop"
  )


form <- norm ~ (a / (1 + exp((b - time) / c))) + (d / (1 + exp((e - time) / f)))

df_mod <- df_test %>%
  left_join(df_test_sum) %>%
  group_by(across(all_of(group_list))) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      nls(form, data = .x,
        start = list(
          a = .x$peak_norm[1],
          b = .x$time_to_growth_mid[1],
          c = .x$growth_scale[1],
          d = .x$peak_decay[1],
          e = .x$time_to_decay_mid[1],
          f = .x$decay_scale[1]
        ),
        algorithm = "port",
        lower = c(a = 0, b = 0, c = 0.01, d = -Inf, e = .x$time_to_growth_mid[1], f = 0.1),
        upper = c(a = Inf, b = .x$max_time[1], c = 10, d = Inf, e = .x$max_time[1], f = 300)
      ),
      error = function(e) NULL
    ))
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

df_unmod <- df_mod %>%
  filter(map_lgl(model, is.null))

df_unmod %>% 
  as.data.frame() %>%
  select(-c(data, model)) %>%
  skim()

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

df_results %>%
  ungroup() %>%
  na.omit() %>%
  filter(time < 68) %>%
  select(sample, assay, time, norm, pred, growth, decay, resid) %>%
  pivot_longer(c(pred, growth, decay), names_to = "series") %>%
  summarize(
    across(c(value, norm, resid), list(mean = mean, sd = sd, min = min, max = max)),
    n = n(),
    .by = c(time, sample, assay, series, )
  ) %>%
  arrange(time) %>%
  ggplot(aes(time, value_mean, color = series)) +
  geom_point(aes(y=norm_mean), size=0.1, color="black") +
  geom_point(aes(y=resid_mean), size=0.1, color="black") +
  geom_line(linewidth = 1, linetype="dashed") +
  scale_color_manual(values = c("darkgreen", "red", "blue")) +
  facet_grid(vars(sample), vars(assay)) +
  labs(x = "Time", y = "Normalized RFU")

# Plot histogram of residuals
df_results %>%
  ggplot(aes(resid)) +
  geom_histogram() +
  scale_x_continuous(limits = c(-1, 1))

df_results %>%
  ggplot(aes(time, resid)) +
  stat_bin_hex(bins = 200, color = NULL) +
  # geom_point() +
  # geom_smooth(method = "loess", se = TRUE) +
  scale_fill_gradient2(low = "darkblue", mid="darkorange", high="darkred", midpoint = 200) +
  facet_grid(cols=vars(assay)) +
  scale_y_continuous(limits = c(-4, 4)) +
  scale_x_continuous(limits = c(0, 72), breaks = seq(0, 70, 4)) +
  coord_cartesian(expand = FALSE) +
  labs(x = "Time", y = "Residuals")
