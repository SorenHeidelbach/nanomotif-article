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
  "ggnewscale"
)
source("src/constants.R")
source("src/themes.R")
source("src/utility.R")

get_axis_order <- function(wide_df) {
  set.seed(1)
  feature_matrix = as.matrix(wide_df[,-1] )
  y_row_order <- pheatmap::pheatmap(as.matrix(feature_matrix), silent = TRUE)$tree_row$order
  y_factor_order <- wide_df[[1]][y_row_order]
  x_row_order <- pheatmap::pheatmap(t(as.matrix(feature_matrix)), silent = TRUE)$tree_row$order
  x_factor_order <- colnames(wide_df[,-1])[x_row_order]
  return(list(
    y = y_factor_order,
    x = x_factor_order
  ))
}


full_genome_benchmark <- fread("data/benchmarks/monocultures_full_genome/Train/Train/sample_performance.tsv")[
  Tool != "nanomotif_read_level"
][
  , sample := `Sample ID`
][
  , precision := true_positives / (true_positives + false_positives)
][
  , recall := true_positives / (true_positives + false_negatives)
][
  , f1 := 2 * (precision * recall) / (precision + recall)
][
  , Best := f1 == max(f1, na.rm = TRUE), by = `Sample ID` 
][
  , sample := str_remove(sample, "CVM\\d+_") %>% str_replace_all(
    "_", " "
  ) %>% str_replace(
    "e coli k12", "Escherichia coli"
  ) %>% str_replace(
    "m ruber", "Meiothermus ruber"
  ) %>%  str_replace(
    "p thermoglucosidasius", "Parageobacillus thermoglucosidasius"
  )
]


# Sample level barplots
tool_order <- c("nanomotif", "modkit", "microbemod", "motifmaker")
sample_level_bar_plot <- ggplot(full_genome_benchmark) +
  aes(x = factor(Tool, levels = tool_order), y = f1) +
  geom_col(aes(fill = Best), position = position_dodge2()) +
  facet_wrap(~sample) +
  geom_text(aes(label = round(f1, 2)), nudge_y = 0.04) +
  theme_bw() +
  scale_fill_manual(values = c("grey60", "forestgreen")) +
  theme(
    axis.text.y = element_markdown(),
    axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_blank(),   # Remove panel background,
    strip.placement = "outside"
  ) +
  labs(
    x = "Tool",
    y = "F1-score",
    fill = "Highest Score"
  )

ggsave(
  "figures/benchmark1_full-genome.png",
  sample_level_bar_plot,
  width = 7,
  height = 5
)








tool_order <- c("nanomotif", "modkit", "microbemod", "motifmaker")
full_genome_benchmark_summary <- fread("data/benchmarks/monocultures_full_genome/Train/Train/summary_performance.tsv")[
  Tool %in% tool_order
]
full_genome_benchmark_summary_long <- melt(full_genome_benchmark_summary, id.vars=c("Tool"), measure.vars = c("f1_score", "precision", "sensitivity"))
full_genome_benchmark_summary_long[
  , Best := value == max(value, na.rm = TRUE), by = variable
][
  , variable := fcase(
    variable == "f1_score", "F1-score",
    variable == "precision", "Precision",
    variable == "sensitivity", "Sensitivity"
  )
]

summary_bar_plot <- ggplot(full_genome_benchmark_summary_long) +
  aes(x = factor(Tool, levels = tool_order), y = value) +
  facet_grid(~variable)  +
  geom_col(aes(fill = Tool)) +  
  geom_text(aes(label = round(value, 2)), nudge_y = 0.03, size = 3) +
  theme_bw() +
  scale_fill_manual(values = c("#ef9321", "#69b3a2","#fbd563", "#403f81"), guide = NULL) +
  theme(
    axis.text.y = element_markdown(),
    axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_blank(),   # Remove panel background,
    strip.placement = "outside"
  ) +
  labs(
    x = "Tool",
    y = "Score",
    fill = "Highest Score"
  )

ggsave(
  "figures/benchmark1_full-genome_summary.png",
  summary_bar_plot,
  width = 7,
  height = 4
)







