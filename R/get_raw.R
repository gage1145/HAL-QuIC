get_raw <- function(file, ...) {
  library(quicR)
  library(tidyverse)

  get_quic_args <- formals(get_quic)
  params <- list(...)
  
  for (param in names(params)) {
    if (param %in% names(get_quic_args)) {
      assign(param, get_quic_args[[param]])
    }
  }

  file %>%
    get_quic(
      transpose_table = get_quic_args$transpose_table,
      norm_point = get_quic_args$norm_point,
      which_table = get_quic_args$which_table,
      window_size = get_quic_args$window_size,
      smooth = get_quic_args$smooth,
      smooth_factor = get_quic_args$smooth_factor,
      zero = get_quic_args$zero,
      by = get_quic_args$by,
      plate = get_quic_args$plate,
      sheet = get_quic_args$sheet,
    )
}
