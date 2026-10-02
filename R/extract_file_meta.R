extract_file_meta <- function(file) {  
  delim_count <- str_count(file, "/")
  rxn <- str_split_i(file, "/", delim_count + 1) %>%
    str_remove("\\.[[:lower:]]+$")
  print(rxn)
  rxn_split <- str_split_1(rxn, "_")
  meta_vec <- c(
    rxn = rxn,
    date = parse_date_time(rxn_split[1], "%Y%m%d"),
    reader = rxn_split[2],
    tech = rxn_split[3]
  )
}
