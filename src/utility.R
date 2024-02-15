

extract_species_from_gtdbtk <- function(str) {
  split_str <- strsplit(str, ";")[[1]]
  species <- strsplit(split_str[length(split_str)], "__")[[1]][2]
  return(species)
}


regex_to_iupac <- function(regex) {
  # Define the mapping from nucleotide combinations to IUPAC codes
  mapping <- c(A = "A", T = "T", C = "C", G = "G",
               AG = "R", GA = "R", CT = "Y", CG = "S", AT = "W",
               GT = "K", AC = "M", CGT = "B", AGT = "D",
               ACT = "H", ACG = "V")
  
  # Function to find and replace bracketed sequences with their IUPAC codes
  process_matches <- function(regex) {
    matches <- gregexpr("\\[([ATCG]+)\\]", regex)[[1]]
    match.lengths <- attr(matches, "match.length")
    
    # If no matches found, return the regex as is
    if(matches[1] == -1) {
      return(regex)
    }
    
    # Process each match found
    for (i in length(matches):1) {
      start <- matches[i]
      end <- start + match.lengths[i] - 1
      match_str <- substr(regex, start + 1, end - 1) # Exclude brackets
      sorted_match_str <- paste(sort(unlist(strsplit(match_str, ""))), collapse = "")
      iupac_code <- mapping[sorted_match_str]
      
      # Replace the match in the original string
      regex <- paste0(substr(regex, 1, start - 1), iupac_code, substr(regex, end + 1, nchar(regex)))
    }
    
    return(regex)
  }
  
  # Process bracketed sequences
  processed_seq <- process_matches(regex)
  
  # Replace dots with 'N'
  final_seq <- gsub("\\.", "N", processed_seq)
  
  return(final_seq)
}




reverse_complement_iupac <- function(dna_sequence) {
  if (is.na(dna_sequence) | is.null(dna_sequence)){
    return(NA)
  }
  # Define the complement for each IUPAC nucleotide symbol
  complement_map <- list(A = "T", T = "A", C = "G", G = "C",
                         R = "Y", Y = "R", S = "S", W = "W",
                         K = "M", M = "K", B = "V", V = "B",
                         D = "H", H = "D", N = "N")
  
  # Reverse the DNA sequence
  reversed_seq <- rev(strsplit(dna_sequence, "")[[1]])
  
  # Replace each nucleotide with its complement according to the IUPAC codes
  complement_seq <- sapply(reversed_seq, function(nucleotide) {
    complement_map[[nucleotide]]
  })
  
  # Collapse the complement sequence back into a single string
  complement_str <- paste(complement_seq, collapse = "")
  
  return(complement_str)
}












extraxt_gtdb_tax <- function(tax_strings) {
  # Initialize lists to store the results for each taxonomic level
  result_list <- vector("list", length = 7)
  names(result_list) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  
  # Loop through each taxonomic string
  for (i in seq_along(tax_strings)) {
    tax_string <- tax_strings[i]
    if (is.na(tax_string)) {
      # Assign NA to all taxonomic levels if the input string is NA
      for (level in names(result_list)) {
        result_list[[level]][i] <- NA
      }
    } else {
      # Split the string into parts based on the semicolon separator
      parts <- strsplit(tax_string, ";")[[1]]
      
      # Initialize a temporary list to store found taxonomic levels
      temp_list <- setNames(vector("list", length(parts)), rep(NA, length(parts)))
      
      for (part in parts) {
        split_part <- strsplit(part, "__")[[1]]
        code <- split_part[1]
        value <- split_part[2]
        
        # Determine the taxonomic level based on the code and assign the value
        level_name <- switch(code,
                             "d" = "domain",
                             "p" = "phylum",
                             "c" = "class",
                             "o" = "order",
                             "f" = "family",
                             "g" = "genus",
                             "s" = "species",
                             NA)
        temp_list[[level_name]] <- value
      }
      
      # Fill in the result list with the values from temp_list or NA if not present
      for (level in names(result_list)) {
        result_list[[level]][i] <- ifelse(!is.null(temp_list[[level]]), temp_list[[level]], NA)
      }
    }
  }
  
  # Convert the list of lists into a list of vectors
  result_vectors <- lapply(result_list, unlist)
  
  return(result_vectors)
}

