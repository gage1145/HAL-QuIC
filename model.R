library(tidyverse)
library(arrow)
library(magrittr)
library(modelr)
library(quicR)
library(zoo)



raw_file <- "data/data.parquet"

group_list <- c("sample", "wells", "dilutions", "assay", "reaction")

df_ <- raw_file %>%
  read_parquet() %>%
  mutate(
    across(c(reaction, wells, dilutions, assay, mortem, sample_type, animal, sample), as.factor)
  ) %>%
  select(-c(norm, deriv))

df_p <- df_ %>%
  filter(sample == "P" & assay == "RT-QuIC" & reaction == levels(reaction)[1]) 
  # group_by(across(all_of(group_list))) %>%
  # mutate(norm = rollmean(norm, 10, na.pad=TRUE)) %>%
  # ungroup()

# df_ttm <- df_p %>%
#   filter(norm == max(norm, na.rm=TRUE), .by = group_list) %>%
#   rename(time_to_max = time)

# df_calcs <- df_p %>%
#   calculate_metrics(
#     group_list,
#     threshold=4, time_col="time", ttt_values="norm", auc_values = "norm", norm_col = "norm", deriv_col = "deriv"
#   ) %>%
#   left_join(df_ttm)


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
    time_to_growth_max   = time[which(norm == peak_norm)],
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

df_test %>%
  left_join(df_test_sum) %>%
  group_by(across(all_of(group_list))) %>%
  mutate(
    growth = peak_norm / (1 + exp((time_to_growth_mid - time) / growth_scale)),
    decay  = peak_decay / (1 + exp((time_to_decay_mid - time) / decay_scale)),
    combined = growth + decay
  ) %>%
  select(wells, time, norm, growth, decay, combined) %>%
  pivot_longer(c(norm, growth, decay, combined), names_to="data") %>%
  mutate(data = factor(data, levels=c("norm", "growth", "decay", "combined"))) %>%
  arrange(data) %>%
  ggplot(aes(time, value, color=data)) +
  geom_line() +
  scale_color_manual(values=c("black", "red", "blue", "green")) +
  # scale_linewidth_manual(values=c(1, 1, 1, 3)) +
  # scale_linetype_manual(values = c(NULL, "dashed", "dashed", NULL)) +
  # geom_line(aes(y=growth), color="red", linetype="dashed") +
  # geom_line(aes(y=decay), color="blue", linetype="dashed") +
  # geom_line(aes(y=combined), color="green", linewidth=1) +
  facet_grid(vars(wells)) +
  scale_x_continuous(breaks=seq(0, 72, 2))



# Fits a double-sigmoid model (growth + decay) to one group's smoothed time series.
# Returns df with added columns: pred, growth, decay, comb.
fit_double_sigmoid <- function(df) {

  max_time    <- max(df$time, na.rm=TRUE)
  max_val     <- df$norm[which.min(df$deriv2)]
  time_to_max <- df$time[which(df$norm == max_val)[1]]
  time_to_mid <- df$time[which.max(df$deriv)]

  df_growth <- df %>%
    # filter(time <= time_to_max)
    mutate(norm = ifelse(time > time_to_max, max_val, norm))

  df_decay <- df %>%
    # filter(time >= time_to_max) %>%
    mutate(
      norm = ifelse(time <= time_to_max, 0, norm - max_val)
      # norm = norm - max_val
    )

  equillibrium     <- min(df_decay$norm, na.rm = TRUE)
  decay_mid <- (max_time - time_to_max) / 2

  form <- norm ~ SSlogis(time, Asym, xmid, scal)

  growth_coefs <- decay_coefs <- cc <- c()

  tryCatch(
    {
      growth_coefs <- nls(
        form, df_growth,
        start = list(Asym = max_val, xmid = time_to_mid, scal = time_to_mid / 2)
      ) |> coef()
      a <- growth_coefs[1]; b <- growth_coefs[2]; c <- growth_coefs[3]
    }, 
    error = function(e) warning("Growth model failed: ", conditionMessage(e))
  )

  tryCatch(
    {
      decay_coefs  <- nls(
        form, df_decay,
        start = list(Asym = equillibrium, xmid = decay_mid, scal = decay_mid / 2)
      ) |> coef()
      d <- decay_coefs[1];  e <- decay_coefs[2];  f <- decay_coefs[3]
    }, 
    error = function(e) warning("Decay model failed: ", conditionMessage(e))
  )

  if (length(growth_coefs) == 0 | length(decay_coefs) == 0) return(
    mutate(df, pred = NA_real_, growth = NA_real_, decay = NA_real_, comb = NA_real_)
  )
  
  tryCatch(
    {
      comb_mod <- nls(
        norm ~ (a / (1 + exp((b - time) / c))) + (d / (1 + exp((e - time) / f))),
        df,
        start = list(a=a, b=b, c=c, d=d, e=e, f=f)
      )
      cc <- coef(comb_mod)
    },
    error = function(e) warning("Combined model failed: ", conditionMessage(e))
  )

  if (length(cc) == 0) return(
    mutate(df, pred = NA_real_, growth = NA_real_, decay = NA_real_, comb = NA_real_)
  )

  tryCatch(
    {
      df %>%
        add_predictions(comb_mod) %>%
        mutate(
          growth = cc[1] / (1 + exp((cc[2] - time) / cc[3])),
          decay  = cc[4] / (1 + exp((cc[5] - time) / cc[6])),
          comb   = growth + decay
        ) %>%
        return()
    }, 
    error = function(e) {
      warning("Failed to add predictions: ", conditionMessage(e))
      mutate(df, pred = NA_real_, growth = NA_real_, decay = NA_real_, comb = NA_real_)
    }
  )
}

df_modeled <- df_test %>%
  group_by(across(all_of(group_list))) %>%
  group_modify(~ fit_double_sigmoid(.)) %>%
  ungroup()

df_unmodeled <- df_modeled %>%
  filter(is.na(pred))


# --- Plots ---

# Equation
# norm ~ (a / (1 + exp((b - time) / c))) + (d / (1 + exp((e - time) / f)))

df_unmodeled %>%
  # filter(wells == "A01") %>%
  ggplot(aes(time)) +
  geom_point(aes(y = norm)) +
  geom_line(aes(y = pred),   color = "red",    linewidth = 1) 
  # geom_line(aes(y = growth), color = "blue",   linewidth = 1) +
  # geom_line(aes(y = decay),  color = "orange", linewidth = 1)

df_modeled %>%
  mutate(resid = norm - pred) %>%
  ggplot(aes(resid)) +
  geom_histogram()
