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
    across(c(reaction, wells, dilutions, assay, mortem, sample_type, animal, sample), as.factor),
    norm = norm - 1
  )

df_p <- df_ %>%
  filter(sample == "P") %>%
  group_by(sample, wells, dilutions, assay, reaction) %>%
  mutate(norm_smooth = rollmean(norm, 10, na.pad=TRUE)) %>%
  ungroup()
#   mutate(norm_smooth = smooth.spline(time, norm, spar=0.4)) %>%
#   add_predictions(loess(norm ~ time, ., span=0.09), var="norm_smooth")

df_ttm <- df_p %>%
  filter(norm_smooth == max(norm_smooth, na.rm=TRUE), .by = group_list) %>%
  rename(time_to_max = time)

df_calcs <- df_p %>%
  calculate_metrics(
    group_list,
    threshold=4, time_col="time", ttt_values="norm", auc_values = "norm", norm_col = "norm", deriv_col = "deriv"
  ) %>%
  left_join(df_ttm)


# Fits a double-sigmoid model (growth + decay) to one group's smoothed time series.
# Returns df with added columns: pred, growth, decay, comb.
fit_double_sigmoid <- function(df) {
  form <- norm_smooth ~ SSlogis(time, Asym, xmid, scal)

  time_to_max <- df$time[which.max(df$norm_smooth)]
  max_val     <- max(df$norm_smooth, na.rm = TRUE)

  df_growth <- df %>%
    mutate(norm_smooth = ifelse(time > time_to_max, max_val, norm_smooth))

  df_decay <- df %>%
    mutate(norm_smooth = ifelse(time < time_to_max, 0, norm_smooth - max_val))

  tryCatch({
    growth_coefs <- coef(nls(form, df_growth))
    decay_coefs  <- coef(nls(form, df_decay))

    a <- growth_coefs[1]; b <- growth_coefs[2]; c <- growth_coefs[3]
    d <- decay_coefs[1];  e <- decay_coefs[2];  f <- decay_coefs[3]

    comb_mod <- nls(
      norm_smooth ~ SSlogis(time, a, b, c) + SSlogis(time, d, e, f),
      df,
      start = list(a=a, b=b, c=c, d=d, e=e, f=f)
    )
    cc <- coef(comb_mod)

    df %>%
      add_predictions(comb_mod) %>%
      mutate(
        growth = cc[1] / (1 + exp((cc[2] - time) / cc[3])),
        decay  = cc[4] / (1 + exp((cc[5] - time) / cc[6])),
        comb   = growth + decay
      )
  }, error = function(e) {
    warning("Model failed: ", conditionMessage(e))
    df %>% mutate(pred = NA_real_, growth = NA_real_, decay = NA_real_, comb = NA_real_)
  })
}

df_modeled <- df_p %>%
  group_by(across(all_of(group_list))) %>%
  group_modify(~ fit_double_sigmoid(.x)) %>%
  ungroup()


# --- Plots ---

# Equation
# norm_smooth ~ (a / (1 + exp((b - time) / c))) + (d / (1 + exp((e - time) / f)))

df_modeled %>%
  # filter(wells == "A01", reaction == levels(reaction)[1]) %>%
  ggplot(aes(time)) +
  geom_point(aes(y = norm_smooth)) +
  geom_line(aes(y = pred),   color = "red",    linewidth = 1) +
  geom_line(aes(y = growth), color = "blue",   linewidth = 1) +
  geom_line(aes(y = decay),  color = "orange", linewidth = 1)

df_modeled %>%
  mutate(resid = norm_smooth - pred) %>%
  ggplot(aes(resid)) +
  geom_histogram()
