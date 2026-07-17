fit_model <- function(
  data, peak_norm, time_to_growth_mid, growth_scale, peak_decay, time_to_decay_mid, decay_scale, max_time, ..., 
  single_only = FALSE, algorithm = "port", weighted = FALSE, irls_iter = 3, irls_s0 = 0.072, irls_k = 0.0094, irls_theta = 1
) {
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

  peak_scalar <- 3
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