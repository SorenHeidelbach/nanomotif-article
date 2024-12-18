```{r}
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
  "openxlsx"
)
source("src/constants.R")
source("src/themes.R")
source("src/utility.R")
```


```{r}

load_bin_info <- function(path, sample){
  assembly_stats <- fread(path)[, c("bin", "Completeness", "Contamination")][, Sample := sample]
  return(assembly_stats)
}
bin_info <- rbind(
  load_bin_info("data/anaerobic_digestor/mmlong2_lite/results/mmlong2_lite_bins.tsv", "anaerobic_digestor"),
  load_bin_info("data/soil_mfd02199/mmlong2_lite/results/mmlong2_lite_bins.tsv", "soil"),
  load_bin_info("data/fecal/mmlong2_lite/results/mmlong2_lite_bins.tsv", "fecal_complex"),
  load_bin_info("data/PaPr00000216MP_nm_0.1.19/NP-PaPr00000216MP_bins.tsv", "fecal_simple")
)

load_bin_consensus <- function(path, sample){
  bin_consensus <- fread(path)[, Sample := sample]
  return(bin_consensus)
}
bin_consensus_motifs <- rbind(
  load_bin_consensus("data/anaerobic_digestor/nanomotif/bin-motifs.tsv", "anaerobic_digestor"),
  load_bin_consensus("data/soil_mfd02199/nanomotif/bin-motifs.tsv", "soil"),
  load_bin_consensus("data/fecal/nanomotif/bin-motifs.tsv", "fecal_complex"),
  load_bin_consensus("data/PaPr00000216MP_nm_0.1.19/nanomotif/bin-motifs.tsv", "fecal_simple")
)

load_mtase_links <- function(path, sample){
  bin_consensus <- fread(path)[, Sample := sample]
  return(bin_consensus)
}

bin_mtase_links <- rbind(
  load_bin_consensus("data/anaerobic_digestor/mtase_linker/nanomotif_assignment_table.tsv", "anaerobic_digestor")),
  load_bin_consensus("data/soil_mfd02199/mtase_linker/nanomotif_assignment_table.tsv", "soil"),
  load_bin_consensus("data/fecal/mtase_linker/nanomotif_assignment_table.tsv", "fecal_complex"),
  load_bin_consensus("data/PaPr00000216MP/mtase_linker/nanomotif_assignment_table.tsv", "fecal_simple"),
)

gathered <- merge(bin_consensus_motifs, bin_mtase_links, all=TRUE) %>% 
  merge(bin_info[, c("Sample", "bin", "Contamination", "Completeness")], all=TRUE, by = c("Sample", "bin"))
gathered[, has_motif := !is.na(motif)]
setorderv(gathered, cols = c("has_motif", "bin", "linked"), order =  -1)
setnames(
  gathered, 
  c("bin", "has_motif", "linked"),
  c("Bin", "Has_motif", "Linked_MTase")
)
gathered <- gathered[, c(1, 2, 17, 16, 15, 4, 5, 3, 6:14)]


write.xlsx(split(gathered[, !"Sample"], gathered$Sample), 'bin_consensus_motifs.xlsx')

