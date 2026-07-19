estimate_params <- function(df) {
  df %>%
    mutate(
      max_time             = map_dbl(data, \(x) max(x$time, na.rm=TRUE)),
      growth_scale         = map_dbl(data, \(x) max(x$deriv, na.rm=TRUE)),
      growth_scale         = map_dbl(data, \(x) max(x$deriv, na.rm=TRUE)),
      peak_norm            = map_dbl(data, \(x) x$norm[which.min(x$deriv2)]),
      time_to_growth_max   = map2_dbl(data, peak_norm, \(x, p) x$time[x$norm == p][1]),
      time_to_growth_mid   = map_dbl(data, \(x) x$time[which.max(x$deriv)]),
      max_equillibrium     = map2_dbl(data, time_to_growth_max, \(x, ttgm) max(x$norm[x$time > ttgm], na.rm=TRUE)),
      min_equillibrium     = map2_dbl(data, time_to_growth_max, \(x, ttgm) min(x$norm[x$time > ttgm], na.rm=TRUE)),
      max_decay            = max_equillibrium - peak_norm,
      min_decay            = min_equillibrium - peak_norm,
      equillibrium         = ifelse(abs(min_decay) > max_decay, min_equillibrium, max_equillibrium),
      time_to_equillibrium = map2_dbl(data, equillibrium, \(x, e) x$time[x$norm == e][1]),
      peak_decay           = equillibrium - peak_norm,
      time_to_decay        = time_to_equillibrium - time_to_growth_max,
      time_to_decay_mid    = time_to_growth_max + time_to_decay / 2,
      decay_slope          = replace_na(peak_decay / time_to_decay, 0),
      decay_scale          = abs(decay_slope)
    )
}