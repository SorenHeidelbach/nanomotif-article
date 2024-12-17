

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


######################## BINNARY #################################
add_qc_to_bin_contamination <- function(contamination, qc, contig_length, sample_name="mmlong2_lite.") {
  df <- qc %>% 
    select(bin:Contamination, Total_Contigs) %>% 
    right_join(contamination) %>% 
    mutate(
      mag_quality = case_when(
        completeness_checkm2 >= 90 & contamination_checkm2 <= 5 ~ "HQ",
        completeness_checkm2 >= 50 & contamination_checkm2 <= 10 ~"MQ"
      ),
      mag_quality = factor(mag_quality, levels = c("HQ", "MQ"))
    ) %>% 
    arrange(mag_quality, Total_Contigs) %>% 
    mutate(
      bin = str_remove(bin, sample_name),
      bin_contig_compare = str_remove(bin_contig_compare, sample_name)
    ) %>% 
    left_join(contig_length)
  
  n_contamination_found_in_bin <- df %>% 
    group_by(bin) %>% 
    summarise(
      n_contaminants = n()
    )
  
  df <- left_join(df, n_contamination_found_in_bin)
  return(df)
} 

filter_comparison_df <- function(comparison, contig_length, sample_name = "mmlong2_lite.")  {
  df <- comparison %>% 
    mutate(
      bin = str_remove(bin, sample_name),
      bin_contig = str_remove(bin_contig, sample_name),
      contig_bin = str_remove(contig_bin, sample_name),
      # Extract motif before the underscore
      motif = str_extract(motif_mod, "^[^_]+"),
      # Extract modification type after the underscore and before the dash
      mod_type = str_extract(motif_mod, "(?<=_)[^-]+"),
      # Extract modification position after the dash, converting to numeric
      mod_position = as.numeric(str_extract(motif_mod, "(?<=-)[0-9]+")),
      motif_axis = paste0(
        str_sub(motif, 1, mod_position),
        "<strong>",
        str_sub(motif, mod_position+1, mod_position+1),
        "<sub><sup>",
        map(mod_type, function(x) MOD_TYPE_PRETTY[[x]]) %>%  unlist(),
        "</sup></sub>",
        "</strong>",
        str_sub(motif, mod_position+2, -1)
      )
    ) %>% 
    left_join(contig_length)
  
  return(df)
}

prepare_bin_consensus <- function(consensus, sample_name = "mmlong2_lite.") {
  df <- consensus %>% 
    mutate(bin = str_remove(bin, sample_name)) %>% 
    mutate(motif_mod = paste0(motif, "_", mod_type,"-", mod_position))
  
  return(df)
}

calculateGCContent <- function(filePath) {
  sequences <- read.fasta(file = filePath) # Read the assembly file
  gcContent <- sapply(sequences, function(seq) {
    GC(seq) # Use the GC function to calculate GC content directly
  })
  
  # Create a data frame with contig names and their GC content
  gcDataFrame <- data.frame(
    contig = names(gcContent),
    gc = gcContent
  )
  
  return(gcDataFrame) # Return the data frame
}

# Check if the sequence has n stretches of m consecutive characters
has_n_character_stretches_of_length_m <- function(sequence, n, m, character) {
  regex_str <- paste0("(", character, "){", m, ",}")
  matches <- str_extract_all(sequence, regex_str)[[1]]
  return(length(matches) >= n)
}

# Function to determine motif type
motif_type <- function(motif_str) {
  # Check for ambiguous motifs
  if (has_n_character_stretches_of_length_m(motif_str, 2, 2, "N")) {
    return("ambiguous")
  }
  # Check for bipartite motifs
  else if (str_detect(motif_str, "(N){3,}")) {
    return("bipartite")
  }
  # Check if motif is a palindrome
  else if (reverse_complement(motif_str) == motif_str) {
    return("palindrome")
  }
  # Default to non-palindrome
  else {
    return("non-palindrome")
  }
}

reverse_complement <- function(seq) {
  # Define the complement map including IUPAC nucleotide codes
  complement <- c(
    "A" = "T", "T" = "A", "C" = "G", "G" = "C",
    "R" = "Y", "Y" = "R", "S" = "S", "W" = "W",
    "K" = "M", "M" = "K", "B" = "V", "V" = "B",
    "D" = "H", "H" = "D", "N" = "N"
  )
  
  # Split the sequence into individual characters
  seq_split <- str_split(seq, "")[[1]]
  
  # Replace each base with its complement and reverse the sequence
  rev_comp <- rev(sapply(seq_split, function(base) {
    if (!is.null(complement[[base]])) {
      return(complement[[base]])
    } else {
      stop(paste("Invalid nucleotide:", base))
    }
  }))
  
  # Combine the reversed complement sequence back into a string
  return(paste(rev_comp, collapse = ""))
}

# Generate all possible k-mers of a given alphabet and size
generate_kmers <- function(k, alphabet = c("A", "T", "C", "G")) {
  kmers <- unlist(lapply(1:k, function(i) {
    apply(expand.grid(rep(list(alphabet), i)), 1, paste, collapse = "")
  }))
  return(kmers)
}

# Find the starting indices of a subsequence in a sequence
subseq_indices <- function(subseq, seq) {
  matches <- str_locate_all(seq, subseq)[[1]]
  return(matches[, 1])  # Return the start positions
}


library(Biostrings)

# Function to load FASTA and return sequence lengths
get_fasta_sequence_lengths <- function(fasta_file) {
  # Read the FASTA file
  fasta_data <- readDNAStringSet(fasta_file)
  
  # Extract sequence IDs
  sequence_ids <- names(fasta_data)
  
  # Calculate sequence lengths
  sequence_lengths <- width(fasta_data)
  
  # Create a data frame
  result <- data.frame(
    contig = sequence_ids,
    length = sequence_lengths,
    stringsAsFactors = FALSE
  )
  
  return(result)
}



reverse_complement_func <- function(sequences) {
  # Define the complement for each nucleotide including IUPAC codes
  complement <- c(
    A = "T", T = "A", G = "C", C = "G",
    R = "Y", Y = "R", S = "S", W = "W",
    K = "M", M = "K", B = "V", D = "H",
    H = "D", V = "B", N = "N",
    a = "t", t = "a", g = "c", c = "g",
    r = "y", y = "r", s = "s", w = "w",
    k = "m", m = "k", b = "v", d = "h",
    h = "d", v = "b", n = "n"
  )
  
  # Apply reverse complement to each sequence in the vector
  reverse_complement <- sapply(sequences, function(sequence) {
    # Split the sequence into individual characters
    sequence_chars <- unlist(strsplit(sequence, split = ""))
    
    # Replace each nucleotide with its complement
    complement_chars <- complement[sequence_chars]
    
    # Reverse the complemented sequence
    paste0(rev(complement_chars), collapse = "")
  })
  
  return(reverse_complement)
}
