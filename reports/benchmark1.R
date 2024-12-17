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

MOD_CODE_TO_PRETTY <- c(
  "m" = "5mC",
  "a" = "6mA",
  "21839" = "4mC"
)
x_order = c("nanomotif", "modkit", "microbemod", "motifmaker")
bench <- fread("data/benchmarks/monocultures_full_genome/Benchmark/Benchmark/summary_performance.tsv")
f1 <- ggplot(bench[Tool %in% x_order]) +
  aes(x=factor(Tool, levels = x_order), y = f1_score, fill = Tool) +
  geom_text(aes(label = round(f1_score, 2)), nudge_y = 0.03, size = 3) +
  geom_col() +
  theme_bw() +
  ylab("F1-score") +
  xlab("") +
  scale_fill_manual(values = c("#ef9321", "#69b3a2","#fbd563", "#403f81"), guide = NULL) +
  ylim(0,1)+ ggtitle("Test dataset") +
  labs(
    subtitle = "Full genome motif identification"
  )
sens <- ggplot(bench[Tool %in% x_order]) +
  aes(x=factor(Tool, levels = x_order), y = sensitivity, fill = Tool) +
  geom_text(aes(label = round(sensitivity, 2)), nudge_y = 0.03, size = 3) +
  geom_col() +
  theme_bw()+
  ylab("Recall") +
  scale_fill_manual(values = c("#ef9321", "#69b3a2","#fbd563", "#403f81"), guide = NULL) +
  xlab("")+
  ylim(0,1)
prec <- ggplot(bench[Tool %in% x_order]) +
  aes(x=factor(Tool, levels = x_order), y = precision, fill = Tool) +
  geom_text(aes(label = round(precision, 2)), nudge_y = 0.03, size = 3) +
  geom_col() +
  theme_bw()+
  ylab("Precision") +
  scale_fill_manual(values = c("#ef9321", "#69b3a2","#fbd563", "#403f81"), guide = NULL) +
  xlab("")+
  ylim(0,1)

(f1 + sens + prec + plot_layout(guides = "collect") &  theme(
    axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)
))  %>% 
ggsave(
  "figures/benchmark1_full-genome_summary-Benchmark-data.png",
  .
)

