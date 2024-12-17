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
  "Biostrings",
  "patchwork"
)
getwd()
source("../src/constants.R")
source("../src/themes.R")
source("../src/utility.R")


load_genomad <- function(path) {
  genomad <- fread(path)
  genomad <- melt(genomad, id.vars = "seq_name")[
      , score := max(value), by = seq_name
    ][
      score == value
    ][
      , value := NULL
    ]
  genomad <- unique(genomad, on = "seq_name")
  return(genomad)
}
load_contig_association <- function(path) {
  return(fread(path))
}

count_non_chromosomal_contigs <- function(genomad, association) {
  non_chromosomal_contig = genomad[variable != "chromosome_score"]$seq_name
  return(sum(association$contig %in% non_chromosomal_contig))
}





count_non_chromosomal_contigs(
  load_genomad("/home/bio.aau.dk/lx38ll/dark-science/motif-identification/analysis/fecal2/genomad_mmlong2/asm_pol_lenfilt_aggregated_classification/asm_pol_lenfilt_aggregated_classification.tsv"),
  load_contig_association("/home/bio.aau.dk/lx38ll/dark-science/nanomotif-article/data/fecal/binnary/include_contigs.tsv")
)

count_non_chromosomal_contigs(
  load_genomad("/home/bio.aau.dk/lx38ll/dark-science/motif-identification/analysis/AD/genomad_mmlong2/asm_pol_lenfilt_aggregated_classification/asm_pol_lenfilt_aggregated_classification.tsv"),
  load_contig_association("/home/bio.aau.dk/lx38ll/dark-science/nanomotif-article/data/anaerobic_digestor/binnary/include_contigs.tsv")
)

count_non_chromosomal_contigs(
  load_genomad("/home/bio.aau.dk/lx38ll/dark-science/nanomotif-article/data/PaPr00000216MP_nm_0.1.19/genomad/NP-PaPr00000216MP_assembly_aggregated_classification/NP-PaPr00000216MP_assembly_aggregated_classification.tsv"),
  load_contig_association("/home/bio.aau.dk/lx38ll/dark-science/nanomotif-article/data/PaPr00000216MP_nm_0.1.19/binnary/include_contigs.tsv")
)

