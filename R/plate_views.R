library(quicR)
source("R/get_raw.R")
source("R/extract_file_meta.R")


existing <- list.files("figures/plate_views", full.names = FALSE)
raw_files <- list.files("raw", pattern = "*.xlsx", full.names = TRUE, recursive = TRUE)
raw_file_names <- map(raw_files, extract_file_meta)

save_plate_views <- function(file) {
  file_meta <- extract_file_meta(file)
  save_file <- 
  file %>%
    get_raw() %>%
    plate_view() %>%
    ggtitle(file_meta$rxn)
  ggsave(paste0(file_meta$rxn, ".png"), path="figures/plate_views/", width = 12, height = 8, dpi = 300)
}
