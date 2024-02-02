

extract_species_from_gtdbtk <- function(str) {
  split_str <- strsplit(str, ";")[[1]]
  species <- strsplit(split_str[length(split_str)], "__")[[1]][2]
  return(species)
}