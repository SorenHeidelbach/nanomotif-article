pacman::p_load(
  "data.table",
  "ggplot2",
  "ggforce",
  "stringr",
  "tidyr",
  "dplyr",
  "ggtext",
  "pheatmap",
  "purrr",
  "here",
  "ggnewscale",
  "patchwork",
  "grid",
  "egg",
  "ggh4x"
)
source("src/constants.R")
source("src/themes.R")
source("src/utility.R")

summary <- fread("data/parameter_sweep/full_genome/metrics_summary.tsv")
summary[
  , minimum_kl_divergence := paste0("Min KL-divergence: ", minimum_kl_divergence)
]
f1 <- ggplot(summary) +
  aes(x = threshold_methylation_general, y = min_motif_score, fill = f1_score) +
  geom_tile(col="black") +
  geom_text(aes(label = f1_score), size = 3) +
  scale_fill_gradientn(
    limits = c(0.7, 0.9), 
    colours = c("white",PLOT_COLORS[[4]],  PLOT_COLORS[[5]]), na.value = "white", 
    values = c(0, 0.5, 1)) +
  facet_wrap(~minimum_kl_divergence, ncol = 1) +
  labs(
    x = "General methylation threshold",
    y = "Minimum motif score",
    title = "F1-score",
    fill = ""
  )


prec <- ggplot(summary) +
  aes(x = threshold_methylation_general, y = min_motif_score, fill = precision) +
  geom_tile(col="black") +
  geom_text(aes(label = precision), size = 3) +
  scale_fill_gradientn(
    limits = c(0.7, 0.9), 
    colours = c("white",PLOT_COLORS[[4]],  PLOT_COLORS[[5]]), na.value = "white", 
    values = c(0, 0.5, 1)) +
  facet_wrap(~minimum_kl_divergence, ncol = 1) +
  labs(
    x = "General methylation threshold",
    y = "Minimum motif score",
    title = "Precision",
    fill = ""
  )

recall <- ggplot(summary) +
  aes(x = threshold_methylation_general, y = min_motif_score, fill = recall) +
  geom_tile(col="black") +
  geom_text(aes(label = recall), size = 3) +
  scale_fill_gradientn(
    limits = c(0.7, 0.9), 
    colours = c("white",PLOT_COLORS[[4]],  PLOT_COLORS[[5]]), na.value = "white", 
    values = c(0, 0.5, 1)) +
  facet_wrap(~minimum_kl_divergence, ncol = 1) +
  labs(
    x = "General methylation threshold",
    y = "Minimum motif score",
    title = "Recall",
    fill = ""
  )

parameter_sweep <- f1 + prec + recall + plot_layout(guides = "collect", axes = "collect")


ggsave("figures/parameter_sweep.png",
       width = 6,
       height = 12,
       parameter_sweep)
