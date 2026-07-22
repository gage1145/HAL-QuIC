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
                      algorithm = "port") 
{
  form1 <- norm ~ (S1 / (1 + exp(a1 * (b1 - time))))
  form2 <- as.formula(paste(deparse(form1), "+ (S2 / (1 + exp(a2 * (b2 - time))))"))


  fit_timeout <- 10 # seconds per nls() call
  fit_control <- nls.control(maxiter = 1000, warnOnly = TRUE)
  nls_bounded <- function(...) {
    setTimeLimit(elapsed = fit_timeout, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
    nls(..., control = fit_control, algorithm = algorithm)
  }

  peak_scalar <- 3
  is_negative_decay <- peak_decay < 0

  start_single <- c(
    S1 = peak_norm, 
    a1 = growth_scale, 
    b1 = time_to_growth_mid
  )

  start_double <- start_single |>
    c(
      S2 = peak_decay, 
      a2 = decay_scale,
      b2 = time_to_decay_mid
    )

  lower_single <- c(
    lower_S1 = peak_norm,
    lower_a1 = 0.1,
    lower_b1 = 0
  )

  lower_double <- lower_single |>
    c(
      lower_S2 = ifelse(is_negative_decay, -peak_norm * peak_scalar / 2, 0),
      lower_a2 = 0,
      lower_b2 = -time_to_decay_mid
    )

  upper_single <- c(
    upper_S1 = peak_norm * peak_scalar,
    upper_a1 = 20,
    upper_b1 = max_time
  )

  upper_double <- upper_single |>
    c(
      upper_S2 = ifelse(is_negative_decay, 0, peak_norm * peak_scalar / 2),
      upper_a2 = 5,
      upper_b2 = max_time
    )

  fit_single <- function(...) {
    single_mod <- NULL
    single_mod <- nls_bounded(form1, data = data, start = start_single, lower = lower_single, upper = upper_single) %>%
      try(silent = TRUE)
    return(single_mod)
  }

  if (single_only) {
    return(fit_single())
  }

  fit_double <- function(...) {

    mod <- NULL
    mod <- nls_bounded(form2, data = data, start = start_double, lower = lower_double, upper = upper_double) %>%
      try(silent = TRUE)

    if (is.null(mod)) {
      mod <- try(fit_single())

      if (is.null(mod)) {
        return(NULL)
      }

      coefs <- coef(mod)
      start = list(
        S1 = coefs[1], a1 = coefs[2], b1 = coefs[3],
        S2 = peak_decay, a2 = decay_scale, b2 = time_to_decay_mid
      )

      mod <- nls_bounded(form2, data = data, start = start, lower = lower_double, upper = upper_double) %>%
        try(silent = TRUE)
    }

    return(mod)
  }

  fit_double()
}