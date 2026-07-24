fit_model <- function(
  data, peak_norm, time_to_growth_mid, growth_scale, peak_decay, 
  time_to_decay_mid, decay_scale, max_time, ...,
  single_only = FALSE, algorithm = "port", peak_scalar = 3
) {

  # Set up timeout and control.
  fit_timeout <- 10
  fit_control <- nls.control(maxiter = 1000, warnOnly = TRUE)
  nls_bounded <- function(...) {
    setTimeLimit(elapsed = fit_timeout, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
    nls(..., control = fit_control, algorithm = algorithm)
  }

  # Single sigmoid function.
  fit <- function(form, start, lower, upper, ...) {
    m <- NULL
    m <- nls_bounded(
      form, data = data, start = start, 
      lower = lower, upper = upper
    ) %>%
      try(silent = TRUE)
    return(m)
  }

  form1 <- norm ~ (S1 / (1 + exp(a1 * (b1 - time))))
  form2 <- as.formula(paste(deparse(form1), "+ (S2 / (1 + exp(a2 * (b2 - time))))"))
  is_negative_decay <- peak_decay < 0


# Single sigmoid fit -----------------------------------------------------


  # Initial starting values and parameter bounds.
  start_single <- c(S1 = peak_norm, a1 = growth_scale, b1 = time_to_growth_mid)
  lower_single <- c(S1 = peak_norm, a1 = 0.1, b1 = 0)
  upper_single <- c(S1 = peak_norm * peak_scalar, a1 = 20, b1 = max_time)

  # Avoid double sigmoid fit if only the single is desired.
  if (single_only) {
    return(fit(form1, start_single, lower_single, upper_single))
  }


# Double sigmoid fit -----------------------------------------------------


  # Initial starting values and parameter bounds.
  start_double <- start_single |>
    c(S2 = peak_decay, a2 = decay_scale, b2 = time_to_decay_mid)
  lower_double <- lower_single |>
    c(
      S2 = ifelse(is_negative_decay, -peak_norm * peak_scalar / 2, 0),
      a2 = 0,
      b2 = -time_to_decay_mid
    )
  upper_double <- upper_single |>
    c(
      S2 = ifelse(is_negative_decay, 0, peak_norm * peak_scalar / 2),
      a2 = 5,
      b2 = max_time
    )

  mod <- fit(form2, start_double, lower_double, upper_double)

  if (inherits(mod, "nls")) return(mod)

  mod <- fit(form1, start_single, lower_single, upper_single)

  if (inherits(mod, "nls")) {
    start_double[1:3] <- coef(mod)
    mod <- fit(form2, start_double, lower_double, upper_double)
  }

  return(mod)
}
