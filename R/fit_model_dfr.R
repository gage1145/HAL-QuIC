fit_model_dfr <- function(df, single_only = FALSE) {
  df %>%
  mutate(model = pmap(
    ., fit_model,
    single_only = single_only, .progress = TRUE
  )) %>%
  filter(map_lgl(model, ~ inherits(.x, "nls"))) %>%
  mutate(
    data = map2(model, data, ~ {
      cc <- coef(.x)
      # ci <- tryCatch(
      #   predFit(.x, newdata = .y, interval = "confidence", level = 0.95),
      #   error = function(e) NULL
      # )
      .y %>%
        add_predictions(.x) %>%
        add_residuals(.x) %>%
        mutate(
          # lower = if (is.null(ci)) NA_real_ else ci[, "lwr"],
          # upper = if (is.null(ci)) NA_real_ else ci[, "upr"],
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
    tlagprime = b1 - a1 / 2,
    intercept = (S1 / 2) - a1 * b1,
    # tlag = -intercept / a1,
    tlag = b1 - (S1 / (2 * a1))
  )
}
