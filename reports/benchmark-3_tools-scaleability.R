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
  "gt",
  "ggnewscale"
)
source("src/constants.R")
source("src/themes.R")
source("src/utility.R")

load_snakemake_rule_benchmark <- function(path) {
  benchmark <- fread(path)
  benchmark[
    , path := path
  ]
  return(benchmark)
}

get_nanomotif_benchmark_path <- function(sample_name) {
  return(
    paste0("data/real_communities/", sample_name, "/benchmarks/nanomotif_0.4.16_benchmark.txt")
  )
}

get_modkit_benchmark_path <-  function(sample_name) {
  return(
    paste0("data/real_communities/", sample_name, "/benchmarks/modkit_find_motifs_benchmark.txt")
  )
}
selected_samples <- c(
  "fecal_simple",
  "fecal_inhouse",
  "ZymoFecal",
  "anaerobic_digester"
)


modkit_benchmarks <- selected_samples %>% lapply(
  get_modkit_benchmark_path
) %>% lapply(
  load_snakemake_rule_benchmark
) %>% rbindlist() 
modkit_benchmarks[
  , sample := str_extract(path, "(?<=communities/).*?(?=/benchmark)")
][
  , tool := "modkit"
]


nanomotif_benchmarks <- selected_samples %>% lapply(
    get_nanomotif_benchmark_path
  ) %>% lapply(
    load_snakemake_rule_benchmark
  ) %>% rbindlist() 
nanomotif_benchmarks[
  , sample := str_extract(path, "(?<=communities/).*?(?=/benchmark)")
][
  , tool := "nanomotif"
]


benchmarks <- rbind(
  modkit_benchmarks,
  nanomotif_benchmarks
)[
  , sample := factor(sample, levels = selected_samples)
][
  , peak_uss_mem_gb := max_uss / 1024
][
  , sample := sample %>% str_replace(
    "fecal_simple", "Fecal, simple"
  ) %>% str_replace(
    "fecal_inhouse", "Fecal, complex"
  ) %>%  str_replace(
    "anaerobic_digester", "Anaerobic Digester"
  ) %>% str_replace(
    "mfd02199_backup", "Soil"
  )
][
  , sample := factor(
    sample, levels = c("Fecal, simple", "Fecal, complex", "ZymoFecal", "Anaerobic Digester", "Soil")
  )
][
  , tool := factor(tool, levels = c("nanomotif", "modkit"))
][
  , cpu_hour := cpu_time / (60*60)
][
  , cpu_hour_formatted := signif(cpu_hour, 2)
][
  , relative_time := s / min(s), by = sample
]

benchmarks %>% select("tool", "cpu_hour_formatted", "h:m:s", "sample") %>% group_by(sample) %>% 
  gt(rowname_col = "tool") %>% 
  tab_spanner(
    label = 'Metrics',
    columns = -tool
  ) %>% 
  tab_style(
    style = cell_fill(color = "gray90"),
    locations = cells_row_groups()
  ) %>% 
  cols_label(
    `h:m:s`="Time [h:m:s]",
    cpu_hour_formatted="CPU time [h]",
    tool="Tool"
  ) |> 
  fmt_number(decimals = 1) |> 
  cols_align('right')
# 
