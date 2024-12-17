pacman::p_load(
  "data.table",
  "ggplot2",
  "stringr",
  "tidyr",
  "dplyr",
  "ggtext",
  "pheatmap",
  "purrr",
  "here",
  "ggnewscale",
  "jsonlite",
  "gt",
  "patchwork",
  "Biostrings"
)
source("src/constants.R")
source("src/themes.R")
source("src/utility.R")

fragment_length_levels <- c(
  "10kb",
  "25kb",
  "50kb",
  "100kb",
  "1000kb"
)
coverage_levels <- c(
  "10x",
  "25x",
  "50x",
  "100x"
)
load_metrics_json <- function(path) {
  metrics_json <- path %>% 
    readLines(warn = FALSE) %>% 
    paste(collapse = "") %>% 
    fromJSON()
  metrics = as.data.table(metrics_json$columns$values) %>% 
    select(-contains("data")) %>% 
    select(-contains("bit")) %>% 
    select(-contains("name"))
  setnames(metrics, metrics_json$columns$name)
  
  return(metrics[sample != "CVM74_Cellulophaga_lytica"])
  
}
calculate_metrics <- function(df) {
  df[
    , positives := false_negatives + true_positives
  ][
    , recall := true_positives / positives
  ][
    , predicted_positives := true_positives + false_positives
  ][
    , precision := true_positives / predicted_positives
  ][
    , f1 := round(
      (2*true_positives) / (2*true_positives + false_positives + false_negatives),
      2
    )
  ]
}


read_ground_truth_motifs <- function(ground_truth_motif_info) {
  # Add columns for MethylationTypeCode and ComplementaryMethylationTypeCode
  ground_truth_motif_info <- ground_truth_motif_info[include == "TRUE"]
  ground_truth_motif_info[, MethylationTypeCode := MOD_PRETTY_TO_CODE[MethylationType]]
  ground_truth_motif_info[, ComplementaryMethylationTypeCode := MOD_PRETTY_TO_CODE[ComplementaryMethylationType]]
  
  # Forward motifs
  motifs_fwd <- ground_truth_motif_info[
    !is.na(MethylatedBase) & !is.na(MethylationTypeCode) & (MethylatedBase != "?"), 
    .(
      sample_id = `Sample ID`,
      species = Species,
      strain = Strain,
      motif = RecognitionSequence,
      mod_type = MethylationTypeCode,
      mod_position = as.integer(MethylatedBase) - 1,
      support = support,
      data_role = `Data role`
    )
  ]
  
  # Reverse motifs
  motifs_rev <- ground_truth_motif_info[
    !is.na(ComplementaryMethylatedBase) & !is.na(ComplementaryMethylationTypeCode) & (ComplementaryMethylatedBase != "?"),
    .(
      sample_id = `Sample ID`,
      species = Species,
      strain = Strain,
      motif = reverse_complement_func(RecognitionSequence),
      mod_type = ComplementaryMethylationTypeCode,
      mod_position = nchar(RecognitionSequence) - as.integer(ComplementaryMethylatedBase),
      support = support,
      data_role = `Data role`
    )
  ]  
  combined_motifs <- rbindlist(list(motifs_fwd, motifs_rev), fill = TRUE)
  
  return(combined_motifs)
}

get_sample_level_metrics <- function(path) {
  sample_level_metrics <- load_metrics_json(path)
  
  # Sample level summary
  sample_level_metrics[
    , fragment_size := factor(fragment_size, levels = fragment_length_levels)
  ][
    , coverage := factor(coverage, levels = coverage_levels)
  ]
  calculate_metrics(sample_level_metrics)
  sample_level_metrics[
    , f1_mean := mean(f1),
    by = .(fragment_size, coverage, sample)
  ]
  return(sample_level_metrics)
}

get_motif_level_metrics <- function(path) {
  sample_level_metrics <- get_sample_level_metrics(path)
  motif_level_metrics <- melt(sample_level_metrics, 
                              id.vars = setdiff(names(sample_level_metrics), c("true_positive_motifs", "false_positive_motifs", "false_negative_motifs")), 
                              measure.vars = c("true_positive_motifs", "false_positive_motifs", "false_negative_motifs"),
                              variable.name = "motif_type",
                              value.name = "motif_list")
  
  motif_level_metrics <- motif_level_metrics[, lapply(.SD, as.character), .SDcols = names(motif_level_metrics)]
  
  motif_level_metrics <- motif_level_metrics[, .(motif_mod_pos = unlist(strsplit(as.character(motif_list), ","))), 
                                             by = c(setdiff(names(motif_level_metrics), "motif_list"))]
  
  motif_level_metrics <- motif_level_metrics[
    , motif := map(motif_mod_pos, function(x) str_extract(x, "^[A-Z]*")) %>%  unlist()
  ][
    , mod_type := map(motif_mod_pos, function(x) str_extract(x, "[a-z]")) %>%  unlist()
  ][
    , motif_sequence_type := lapply(motif, function(x) motif_type(x)) %>%  unlist()
  ][
    , .(
      true_positives = as.numeric(sum(motif_type == "true_positive_motifs")),
      false_negatives = as.numeric(sum(motif_type == "false_negative_motifs")), 
      false_positives = as.numeric(sum(motif_type == "false_positive_motifs"))
    ), by = .(fragment_size, coverage, motif_sequence_type)
  ][
    , fragment_size := factor(fragment_size, levels = fragment_length_levels)
  ][
    , coverage := factor(coverage, levels = coverage_levels)
  ]
  motif_level_metrics <- calculate_metrics(motif_level_metrics)
  return(motif_level_metrics)
}

get_summary_metrics <- function(path) {
  sample_level_metrics <- get_sample_level_metrics(path)
  summary_metrics <- sample_level_metrics[
    , .(
      true_positives = sum(true_positives),
      false_negatives = sum(false_negatives), 
      false_positives = sum(false_positives)
    ), by = .(fragment_size, coverage, sample)
  ][
    , positives := false_negatives + true_positives
  ][
    , recall := true_positives / positives
  ][
    , predicted_positives := true_positives + false_positives
  ][
    , precision := true_positives / predicted_positives
  ][
    , f1 := round(
      (2*true_positives) / (2*true_positives + false_positives + false_negatives),
      2
    )
  ][
    , fragment_size := factor(fragment_size, levels = fragment_length_levels)
  ][
    , coverage := factor(coverage, levels = coverage_levels)
  ]
  return(summary_metrics)
}




gold_standard_motifs <- fread("data/gold_standard_motifs.tsv") %>% 
  read_ground_truth_motifs()


motif_expected_count_benchmark <- gold_standard_motifs[data_role == "Benchmark"][
  , motif_sequence_type := lapply(motif, function(x) motif_type(x)) %>%  unlist()
][
  , .(
    count = .N
  ), by = .(motif_sequence_type)
]
motif_expected_count_train <- gold_standard_motifs[data_role == "Train"][
  , motif_sequence_type := lapply(motif, function(x) motif_type(x)) %>%  unlist()
][
  , .(
    count = .N
  ), by = .(motif_sequence_type)
]

gold_standard_motifs <- gold_standard_motifs[data_role == "Benchmark"][
    , fragment_lengths := list(fragment_length_levels)
  ]
gold_standard_motifs <- gold_standard_motifs[, .(fragment_size = unlist(fragment_lengths)), 
                    by = c(setdiff(names(gold_standard_motifs), "fragment_lengths"))]

count_motifs_in_fasta <- function(fasta_file, motif) {
  # Validate input motif using IUPAC symbols
  iupac_symbols <- c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N")
  if (!all(str_split(motif, "")[[1]] %in% iupac_symbols)) {
    stop("The motif contains invalid characters. Use IUPAC nucleotide symbols.")
  }
  
  # Read the FASTA file
  assembly <- readDNAStringSet(fasta_file)
  
  # Translate IUPAC motif to regex pattern
  iupac_to_regex <- c(
    A = "A", C = "C", G = "G", T = "T",
    R = "[AG]", Y = "[CT]", S = "[GC]", W = "[AT]",
    K = "[GT]", M = "[AC]", B = "[CGT]", D = "[AGT]",
    H = "[ACT]", V = "[ACG]", N = "[ACGT]"
  )
  motif_regex <- str_replace_all(motif, setNames(iupac_to_regex, names(iupac_to_regex)))
  
  # Initialize vector to store counts
  counts <- numeric(length(assembly))
  names(counts) <- names(assembly)  # Assign sequence names from the FASTA file
  
  # Count motif occurrences in each sequence and its reverse complement
  for (i in seq_along(assembly)) {
    # Get forward and reverse complement sequences
    forward_seq <- as.character(assembly[i])
    reverse_seq <- as.character(reverseComplement(assembly[i]))
    
    # Count matches on forward strand
    forward_matches <- str_locate_all(forward_seq, motif_regex)[[1]]
    forward_count <- nrow(forward_matches)
    
    # Count matches on reverse complement
    reverse_matches <- str_locate_all(reverse_seq, motif_regex)[[1]]
    reverse_count <- nrow(reverse_matches)
    
    # Store total count for this sequence
    counts[i] <- forward_count + reverse_count
  }
  
  # Return counts for each sequence
  return(counts)
}

gold_standard_motifs[
  , motif_count := mapply(
    function(sample, motif, fragment_size) {
      fasta_path <- paste0("data/monoculture_subset/chunked_assemblies/", fragment_size ,"/", sample, "-", fragment_size, ".fasta")
      count_motifs_in_fasta(fasta_path, motif) %>%  unlist()
    },
    sample = sample_id,
    motif = motif,
    fragment_size = fragment_size
  )
]
gold_standard_motifs_contig <- gold_standard_motifs[
  , .(motif_count = unlist(motif_count))
  , by = c(setdiff(names(gold_standard_motifs), "motif_count"))
][
  , motif_sequence_type := lapply(motif, function(x) motif_type(x)) %>%  unlist()
][
  , fragment_size := factor(fragment_size, levels = fragment_length_levels)
][
  , motif_count_log_friendly := fifelse(
    motif_count == 0, 0.7, motif_count
  )
]

motif_type_count_in_fragments <- gold_standard_motifs_contig[
  , .(
    mean_motif_count = mean(motif_count),
    median = median(motif_count),
    sd = sd(motif_count),
    min = min(motif_count),
    max = max(motif_count),
    q75 = quantile(motif_count, probs = 0.75),
    q25 = quantile(motif_count, probs = 0.25)
  ), by = .(motif_sequence_type, fragment_size)
][
  , plot_motif_sequence_type := fcase(
    motif_sequence_type == "palindrome", "Palindrome<br />e.g. G<strong>A<sub><sup>m6</sup></sub></strong>TC<br />N = 18",
    motif_sequence_type == "non-palindrome", "Non-Palindrome<br />e.g. GGC<strong>C<sub><sup>m4</sup></sub></strong>T<br />N = 16",
    motif_sequence_type == "bipartite", "Bipartite<br />e.g. GG<strong>A<sub><sup>m6</sup></sub></strong>NNNNNCTTA<br />N = 12"
  )
]

ggplot(gold_standard_motifs_contig) +
  aes(x = motif_count_log_friendly, fill = sample_id) +
  geom_density(col="black", alpha = 0.5) +
  ggh4x::facet_grid2(fragment_size~motif_sequence_type, scales = "free", independent = "y") +
  scale_x_log10(n.breaks = 8) +
  geom_vline(
    data=motif_type_count_in_fragments,
    aes(xintercept = median),
    color = "red",
    linewidth = 1
  ) +
  geom_vline(
    data=motif_type_count_in_fragments,
    aes(xintercept = mean_motif_count),
    color = "blue",
    linewidth = 1
  ) +
  theme_bw() 
  


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


dataset <- "Benchmark"
nanomotif_metrics_path <- paste0(
  "/home/bio.aau.dk/lx38ll/dark-science/motif-identification/analysis/monoculture_subset/summary_report_",
  dataset, 
  "/",
  dataset,
  "/nanomotif_chunks_results.json")
modkit_metrics_path <- paste0(
  "/home/bio.aau.dk/lx38ll/dark-science/motif-identification/analysis/monoculture_subset/summary_report_",
  dataset, 
  "/",
  dataset,
  "/modkit_chunks_results.json")


###############################################################################
# Summary of metrics over samples
nanomotif_summary <- get_summary_metrics(nanomotif_metrics_path)[, tool := "Nanomotif"]
modkit_summary <- get_summary_metrics(modkit_metrics_path)[, tool := "Modkit"]

summary_metrics <- rbind(
  nanomotif_summary,
  modkit_summary
)[
  , fragment_size := paste0(fragment_size, "p")
]
summary_plots <- c("f1", "recall", "precision") %>% map(function(x){
  ggplot(summary_metrics) +
  aes(x = fragment_size, y = coverage, fill = get(x)) +
    geom_tile(color = "black")  +
    geom_text(aes(label = round(get(x), 2)), size = 3) +
    facet_grid(tool~sample) +
    scale_fill_gradientn(
      limits = c(0, 1), 
      colours = c("white",PLOT_COLORS[[4]],  PLOT_COLORS[[5]]), na.value = "white", 
      values = c(0, 0.5, 1)) +
    ggtitle(x)

  })
sample_summary_plot <- summary_plots %>%  wrap_plots(ncol = 1)






################################################################################
# Plots split on motif type and modtype
nanomotif_motif_metrics <- get_motif_level_metrics(nanomotif_metrics_path)[, tool := "Nanomotif"]
modkit_motif_metrics <-  get_motif_level_metrics(modkit_metrics_path)[, tool := "Modkit"]


plot_motif_metrics <- function(motif_metrics, title, metric="f1") {
  
  plot_df <- copy(motif_metrics)
  plot_df[
    , plot_motif_sequence_type := fcase(
      motif_sequence_type == "palindrome", "Palindrome<br />e.g. G<strong>A<sub><sup>m6</sup></sub></strong>TC<br />N = 18",
      motif_sequence_type == "non-palindrome", "Non-Palindrome<br />e.g. GGC<strong>C<sub><sup>m4</sup></sub></strong>T<br />N = 16",
      motif_sequence_type == "bipartite", "Bipartite<br />e.g. GG<strong>A<sub><sup>m6</sup></sub></strong>NNNNNCTTA<br />N = 12"
    )
  ]
  plot_df[
    , fragment_size := factor(str_remove(fragment_size, "kb"), levels = str_remove(levels(fragment_size), "kb"))
  ]
  ggplot(plot_df[motif_sequence_type != "ambiguous"]) +
    aes(x = fragment_size, y = coverage, fill = round(get(metric), 2)) +
    geom_tile(color = "black")  +
    geom_text(aes(label = round(get(metric), 2))) +
    facet_wrap(~plot_motif_sequence_type, ncol = 1) +
    scale_fill_gradientn(
      limits = c(0, 1), 
      colours = c("white",PLOT_COLORS[[4]],  PLOT_COLORS[[5]]), na.value = "white", 
      values = c(0, 0.5, 1)) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.grid.major = element_blank(),  # Remove minor gridlines
      panel.background = element_blank(),   # Remove panel background,  
      plot.title = element_text(vjust = 0.5, hjust = 0.5),
      #panel.border = element_rect(colour = "#141414", fill=NA, linewidth = 1),
      strip.text.y.left = element_markdown(angle = 0, vjust = 1, hjust = 0, size =15),
      strip.text.x.top = element_markdown()
      #legend.position="bottom"
    ) +
    labs(
      y = "Coverage",
      x = "Fragment size [kbp]",
      color = "",
      fill = "F1-score",
      title = title
    ) +
    geom_text(
      data = motif_type_count_in_fragments,
      aes(
        x = str_remove(fragment_size, "kb"),
        y = 0,
        label = paste0(median)
      ),
      vjust = 1,
      inherit.aes = FALSE,
      size = unit(3, "pt")
    ) +
    annotate(
      "text",
      y = 0.1,
      x = "50",
      label = "Median motif observations pr. fragment",
      vjust = 0,
      size = unit(3, "pt")
    ) +    
    annotate(
      "text",
      y = -0.5,
      x = "50",
      label = "",
      vjust = 0
    )
}
benchmark2 <- list(
  plot_motif_metrics(nanomotif_motif_metrics, "Nanomotif"),
  plot_motif_metrics(modkit_motif_metrics, "Modkit")
) %>%  wrap_plots()  +
  plot_layout(guides = "collect", axes = "collect")
benchmark2_recall <- list(
  plot_motif_metrics(nanomotif_motif_metrics, "Nanomotif", metric = "recall"),
  plot_motif_metrics(modkit_motif_metrics, "Modkit", metric = "recall")
) %>%  wrap_plots(ncol = 2)  +
  plot_layout(guides = "collect")
benchmark2_precision <- list(
  plot_motif_metrics(nanomotif_motif_metrics, "Nanomotif", metric = "precision"),
  plot_motif_metrics(modkit_motif_metrics, "Modkit", metric = "precision")
) %>%  wrap_plots()  +
  plot_layout(guides = "collect")


benchmark2_recall <- (benchmark2_recall + plot_annotation(
  title = 'Recall',
  theme = theme(plot.title = element_text(size = 32, hjust = 0))
))
benchmark2_precision <- (benchmark2_precision +plot_annotation(
  title = 'Precision',
  theme = theme(plot.title = element_text(size = 32, hjust = 0))
))
ggsave(
  paste0("figures/benchmark2_", dataset, ".png"),
  height = 8,
  width = 7,
  benchmark2
)
ggsave(
  paste0("figures/benchmark2_", dataset, ".pdf"),
  height = 8,
  width = 6.5,
  benchmark2
)
ggsave(
  paste0("figures/benchmark2_", dataset, ".recall.png"),
  height = 8,
  width = 7,
  benchmark2_recall
)
ggsave(
  paste0("figures/benchmark2_", dataset, ".recall.pdf"),
  height = 8,
  width = 7,
  benchmark2_recall
)
ggsave(
  paste0("figures/benchmark2_", dataset, ".precision.png"),
  height = 8,
  width = 7,
  benchmark2_precision
)
ggsave(
  paste0("figures/benchmark2_", dataset, ".precision.pdf"),
  height = 8,
  width = 7,
  benchmark2_precision
)


