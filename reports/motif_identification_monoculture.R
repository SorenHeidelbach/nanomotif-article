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

################################################################################
MOD_PRETTY_TO_CODE <- setNames(names(MOD_CODE_TO_PRETTY), MOD_CODE_TO_PRETTY)

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
      support = support
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
      support = support
    )
  ]  
  combined_motifs <- rbindlist(list(motifs_fwd, motifs_rev), fill = TRUE)
  
  return(combined_motifs)
}

################################################################################
monoculture_data_folder <- "data/monocultures"
monocultures <- list.files(monoculture_data_folder) %>% 
  lapply(function(x) {
    if (!file.exists(file.path(monoculture_data_folder, x))) return(NA)
    else {
      path = file.path(monoculture_data_folder, x, "nanomotif_0.4.16/bin-motifs.tsv")
      if (!file.exists(path)) return(NA)
      fread(path)
    }
  }) %>%  `[`(!is.na(.)) %>% rbindlist()
monocultures[
  , sample_type := "Isolates"
]



################################################################################
mock_data_folder <- "data/mock_communities"
mocks <- list.files(mock_data_folder) %>% 
  lapply(function(x) {
    print(x)
    if (!file.exists(file.path(mock_data_folder, x))) return(NA)
    else {
      path = file.path(mock_data_folder, x, "nanomotif_0.4.16_zymo_ref/bin-motifs.tsv")
      print(path)
      if (!file.exists(path)) return(NA)
      motifs <- fread(path)
      motifs[
        , sample_type := x
      ]
      return(motifs)
    }
  }) %>%  `[`(!is.na(.)) %>% rbindlist()



################################################################################
bin_motifs <- rbind(
  monocultures,
  mocks,
  fill=TRUE
)
bin_motifs_complements <- bin_motifs[
  motif_complement != ""
][
  (motif != motif_complement) | (mod_position != mod_position_complement)
][
  , motif := motif_complement
][
  , mod_position := mod_position_complement
][
  , n_mod_bin := n_mod_complement
][
  , n_nomod_bin := n_nomod_complement
]

bin_motifs <- rbind(
  bin_motifs,
  bin_motifs_complements,
  fill=TRUE
)
bin_motifs[
  , c("motif_complement", "mod_position_complement", "n_mod_complement", "n_nomod_complement") := NULL
]
bin_motifs[
  , identified := TRUE
]



################################################################################
bin_motifs_scored_contigs <- fread(
  "data/motif_scoring/version_0.4.16/bin-motifs-from-mono-and-mock-scored.tsv"
)[
  , group := sample
]

bin_motifs_scored_contigs[
  , group := fcase(
    str_detect(contig, "plasmid"), contig,
    str_detect(contig, "chromosome"), contig,
    sample == "ZymoMock", str_extract(contig, ".*-") %>% str_remove("-"),
    rep(TRUE, length(sample)), sample
  )
][
  sample == "ZymoMock", sample := str_extract(contig, ".*-") %>% str_remove("-")
]

################################################################################
bin_motifs_scored <- bin_motifs_scored_contigs[
  , .(
    n_mod_bin = sum(n_mod),
    n_nomod_bin = sum(n_nomod)
  ), by = .(group, sample, motif, filepath)
][
  , bin := sample
][
  , motif_seq := lapply(motif, function(x) str_split(x, "_")[[1]][[1]]) %>%  unlist()
][
  , mod_type := lapply(motif, function(x) str_split(x, "_")[[1]][[2]]) %>%  unlist()
][
  , mod_position := lapply(motif, function(x) str_split(x, "_")[[1]][[3]] %>%  as.numeric()) %>%  unlist()
][
  , motif := motif_seq
][
  , c("bin", "motif", "mod_type", "mod_position", "n_mod_bin", "n_nomod_bin", "group", "filepath")
][
  mod_type %in% MOD_PRETTY_TO_CODE
]


################################################################################
gold_standard_motifs <- fread("data/gold_standard_motifs.tsv", select = 1:14) %>% 
  read_ground_truth_motifs() 

gold_standard_motifs[
  , gold_standard := TRUE
]


################################################################################
# Merge Directly detected and all motifs scored and gold
bin_motifs_full <- merge(
  bin_motifs,
  bin_motifs_scored,
  by = c("bin", "motif", "mod_type", "mod_position"),
  all = TRUE,
  suffixes = c("", ".scored")
)[
  , bin := group
] %>% merge(
  gold_standard_motifs,
  all = TRUE,
  by.x =c(
    "bin",
    "motif",
    "mod_type",
    "mod_position"
  ),
  by.y = c(
    "sample_id",
    "motif",
    "mod_type",
    "mod_position"
  )
)

bin_motifs_full[
  , n_mod_bin.scored := fifelse(
    is.na(n_mod_bin),
    n_mod_bin.scored,
    n_mod_bin
  )
][
  , n_nomod_bin.scored := fifelse(
    is.na(n_nomod_bin),
    n_nomod_bin.scored,
    n_nomod_bin
  )
][
  , sample_type := fcase(
    str_detect(filepath, "/Zymo[OHG]"), str_extract(filepath, "Zymo[A-Za-z]+"),
    default = "Isolates"
  )
][
  , identified := fifelse(
    is.na(identified),
    "No",
    "Yes"
  )
][
  is.na(gold_standard), gold_standard := FALSE
]

bin_motifs_full[
  , mean := n_mod_bin.scored / (n_mod_bin.scored + n_nomod_bin.scored)
][
  , mod_type_pretty := map(mod_type, function(x) MOD_TYPE_PRETTY[[x]]) %>%  unlist()
][
  , motif_markdown := paste0(
    fifelse(mod_position <= 0, "", str_sub(motif, 1, mod_position)),
    "<strong>",
    str_sub(motif, mod_position+1, mod_position+1),
    "<sub><sup>",
    map(mod_type, function(x) MOD_TYPE_PRETTY[[x]]) %>%  unlist(),
    "</sup></sub>",
    "</strong>",
    str_sub(motif, mod_position+2, -1)
  )
][
  , name := str_remove(bin, "CVM\\d+_") %>% str_replace_all(
    "_", " "
  ) %>% str_remove(
    "DSMZ\\d+-"
  ) %>% str_replace(
    "e coli k12", "Escherichia coli MG1655"
  ) %>% str_replace( 
    "m ruber", "Meiothermus ruber"
  ) %>%  str_replace(
    "p thermoglucosidasius", "Parageobacillus thermoglucosidasius"
  ) %>% str_replace(
    "e coli dcm-dam", "Escherichia coli dam-/dcm-"
  ) %>% str_remove(
    " chromosome"
  ) %>% str_replace(
    "plasmidl", "plasmid1"
  ) %>% str_remove(
    " complete genome"
  ) %>% str_replace(
    "b2207", "B2207"
  ) %>% 
    str_remove(" Concatenation of \\d+ sequences") %>% 
    str_remove("final") %>% 
    str_replace_all("\\.+", " ") %>% 
    str_remove("^[A-Za-z]+ [A-Za-z]+-") %>% 
    str_remove("^[A-Za-z]+ [A-Za-z]+ [A-Z]+\\d+-")
][
  str_detect(name, "Escherichia coli$") & ((sample_type == "ZymoGut") | (sample_type == "ZymoHMW")), name := "Escherichia coli B1109"
][
  str_detect(name, "Escherichia coli plasmid$") & ((sample_type == "ZymoGut") | (sample_type == "ZymoHMW")), name := "Escherichia coli B1109 plasmid"
][
  str_detect(name, "Saccharomyces cerevisiae"), name := "Saccharomyces cerevisiae"
][
  str_detect(name, "coli"), name := str_replace(name, "E.coli", "Escherichia coli") 
][
  str_detect(name, "plasmid"), name := str_replace(name, "plasmid.*", "**\\0**")
][
  , pretty_modtype := fcase(
    mod_type == "a", "6mA",
    mod_type == "m", "5mC",
    mod_type == "21839", "4mC"
  )
][
  , name := str_replace(name, "[A-Za-z]+ [A-Za-z]+", "*\\0*")
][
  , expected_motif := fcase(
    support == "New", "New",
    support != "New", "Expected",
    default = "New"
  )
]

################################################################################
##### Remove motifs with no bin motifs od mean greater than 0.5
bin_motifs_full_filt <- bin_motifs_full[TRUE][
  , keep := any((mean > 0.5) & (identified == "Yes")), by = .(motif, mod_type, mod_position)
][keep == TRUE][, keep := NULL]




################################################################################
# Conig lengths
contig_lengths <- list.files(monoculture_data_folder) %>% 
  lapply(function(x) {
    if (!file.exists(file.path(monoculture_data_folder, x))) return(NA)
    else {
      path = file.path(monoculture_data_folder, x, "assembly.polished.fasta")
      if (!file.exists(path)) return(NA)
      df <- get_fasta_sequence_lengths(path) %>% as.data.table()
      df[, sample := x]
    }
  }) %>%  `[`(!is.na(.)) %>% rbindlist()


################################################################################
# Add read level meth
motifs_scored_contig <- list.files(monoculture_data_folder) %>% 
  lapply(function(x) {
    if (!file.exists(file.path(monoculture_data_folder, x))) return(NA)
    else {
      path = file.path(monoculture_data_folder, x, "pymethylation_util_nanomotif_0.4.16/pymeth_motifs_scored.tsv")
      if (!file.exists(path)) return(NA)
      temp <- fread(path)
      temp[,sample := x]
    }
  }) %>%  `[`(!is.na(.)) %>% rbindlist()
motifs_scored_contig_gs <- list.files(monoculture_data_folder) %>% 
  lapply(function(x) {
    if (!file.exists(file.path(monoculture_data_folder, x))) return(NA)
    else {
      path = file.path(monoculture_data_folder, x, "pymethylation_util_nanomotif_0.4.16/pymeth_motifs_scored_gold_standard.tsv")
      if (!file.exists(path)) return(NA)
      temp <- fread(path)
      temp[,sample := x]
    }
  }) %>%  `[`(!is.na(.)) %>% rbindlist()


motifs_scored_contig <- rbind(motifs_scored_contig, motifs_scored_contig_gs) %>% 
  merge(contig_lengths)

motifs_scored_bin <- motifs_scored_contig[
  , .(
    mean_median = mean(median),
    weighted_median = sum(median * (length/sum(length))),
    max_median = max(median),
    n_contigs = .N,
    mean_mean_read_cov = mean(mean_read_cov),
    n_motif_obs = sum(N_motif_obs)
    ), by = .(sample, motif, mod_type, mod_position)
]





################################################################################
################################################################################
################################################################################
################################################################################

# Seperate plots
isolate_selection_bin <- bin_motifs_full[gold_standard == TRUE]$bin %>%  unique()
isolate_selection <- bin_motifs_full[gold_standard == TRUE]$name %>%  unique()

row1 <- c("*Anabaena variabilis*", "*Kangiella aquimarina*","*Pelobacter carbinolicus*" ,"*Meiothermus ruber*"  ,   "*Salmonella bongori*")
row2 <- c("*Shewanella oneidensis*", "*Desulfobacca acetoxidans*","*Sphaerobacter thermophilus*", "*Escherichia coli* MG1655", "*Cellulophaga lytica*",  "*Thermanaerovibrio acidaminovorans*", "*Zymomonas mobilis*")
plot_data_separated <- bin_motifs_full[name %in% isolate_selection][
  (gold_standard == TRUE) | (identified == "Yes")
][
  !str_detect(name, "natans") 
][
  , gold_standard_chr := fifelse(gold_standard, "yes", "no") 
][
  , rows := fcase(
    name %in% row1, "row1",
    name %in% row2, "row2",
    default = "row3"
  )
] %>%  merge(motifs_scored_bin, by.x = c("bin", "motif", "mod_type", "mod_position"), by.y = c("sample", "motif", "mod_type", "mod_position"), all.x = TRUE)
plot_data_separated[, .N, by = rows]
plot_data_separated[
  , n_motif_obs := fcase(
    nchar(n_motif_obs) < 4, as.character(n_motif_obs),
    nchar(n_motif_obs) < 7, str_replace(n_motif_obs,"...$", "k"),
    nchar(n_motif_obs) < 10, str_replace(n_motif_obs,"......$", "m")
  )
]
plot_data_separated <- unique(plot_data_separated)
plot_row <- function(df) {
  df[
    , motif_sequence_type := lapply(motif, function(x) motif_type(x)) %>%  unlist()
  ]
  setorderv(df, c("expected_motif", "identified", "mod_type", "motif_sequence_type", "mod_position"), order = c(-1, -1, 1, 1, -1))
  df[
    , name := str_replace(name, "Escherichia", "E.") %>% 
      str_remove("MG1655") %>% 
      str_replace("Thermanaerovibrio", "T.") %>% 
      str_remove_all("\\*")
  ]
  df[
    , motif_markdown := fifelse(
      identified == "No",
      paste0("<span style='color:#bc7100'>", motif_markdown, "</span>"),
      paste0("<span style='color:", PLOT_COLORS[[2]], "'>", motif_markdown, "</span>"),
    )
  ]
  tile_border_width = 0
  motif_order <- unique(df$motif_markdown)
  df[
    , support := fifelse(
      is.na(support),
      "Other",
      support
    )
  ][
    , motif_markdown := factor(motif_markdown, levels = motif_order)
  ][
    ,expected_motif := fifelse(expected_motif == "New", "Other", expected_motif)
  ]
  text_size = 2.5
  axis_orders <- get_axis_order(dcast(df, name ~ motif_markdown, value.var = "mean", fill = 0, fun.aggregate = mean))
  p <- ggplot(df) +
    aes(x = motif_markdown, y = 1) +
    
    geom_tile(aes(y=1, fill=round(weighted_median, 2)),  linewidth = tile_border_width, width = 0.9) +
    geom_text(aes(y=1, label = round(weighted_median*100, 0)), size = text_size) +
    scale_fill_gradientn(
      labels = scales::percent,
      guide = "none",
      limits = c(0, 1), 
      colours = c("white",PLOT_COLORS[[4]],  PLOT_COLORS[[5]]), na.value = "white", 
      values = c(0, 0.5, 1)) +
    new_scale("fill") +
    
    geom_tile(aes(y=2, fill=mean),  linewidth = tile_border_width, width = 0.9) +
    geom_text(aes(y=2,  label = round(mean*100, 0)), size = text_size) +
    scale_fill_gradientn(
      labels = scales::percent,
      limits = c(0, 1), 
      colours = c("white",PLOT_COLORS[[4]],  PLOT_COLORS[[5]]), na.value = "white", 
      values = c(0, 0.5, 1)) +
    
    geom_tile(aes(y=3), fill = "white",  linewidth = tile_border_width, width = 0.9) +
    geom_text(aes(y=3, label = n_motif_obs), size = text_size) +
    
    scale_y_continuous(breaks = c(1,2,3), labels=c("Median motif methylation [%]", "Motifs >70% methylated [%]", "Motif genome count")) +
    new_scale("color") +
    geom_text(
      data = df %>% mutate(x = if_else(identified == "Yes", "", " ")),
      aes(col = x, label = x, y = 1), alpha = 0) +
    facet_nested(~name + expected_motif, scale="free_x", space="free_x")  +
    theme_bw() +
    theme(
      axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_markdown(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.grid.major = element_blank(),  # Remove minor gridlines
      panel.background = element_blank(),   # Remove panel background,
      #panel.border = element_rect(colour = "#141414", fill=NA, linewidth = 1),
      strip.text.x = element_text(size =8, face = "italic"),
      strip.text.x.nested = element_text(size =8, face = "italic"),
      legend.position="bottom"
    ) +
    labs(
      y = NULL,
      x = "",
      color = "",
      fill = ""
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          alpha = c(1, 1),
          shape = c(0, 0),
          size = c(5, 5),
          color = c(PLOT_COLORS[[2]], "#bc7100"),
          linetype = c(0, 0),
          linewidth = c(0, 0),
          label = c("Identified", "Not identified")
        )
      )
    )
  return(p)
}
row_plots <-  split(plot_data_separated, by = "rows") %>% lapply(
  plot_row
)

main_figure_plot <- wrap_plots(row_plots, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position="bottom")

ggsave(
  "figures/gold_standard_motif_identification.png",
  width = 15,
  height = 6.5,
  main_figure_plot
)
ggsave(
  "figures/gold_standard_motif_identification.pdf",
  width = 15,
  height = 6.5,
  main_figure_plot
)










################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
######################### Plots not used #######################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

plot_single_isolate <- function(df){
  setorderv(df, c("support", "identified", "mod_type"), order = -1)
  motif_order <- unique(df$motif_markdown)
  df[
    , support := fifelse(
      is.na(support),
      "New",
      support
    )
  ][
    , motif_markdown := factor(motif_markdown, levels = motif_order)
  ]
  p <- ggplot(df) +
    aes(x = motif_markdown) +
    geom_tile(aes(fill = mean, y=1),linewidth = 0.5, width = 0.9, color = "black") +
    geom_text(aes(y=1, label = n_mod_bin.scored)) +
    geom_tile(aes(fill = mean, y=2), linewidth = 0.5, width = 0.9, color = "black") +
    geom_text(aes(y=2, label = n_mod_bin.scored + n_nomod_bin.scored)) +
    geom_tile(aes(fill = mean, y=3), linewidth = 0.5, width = 0.9, color = "black") +
    geom_text(aes(y=3, label = round(mean, 1))) +
    
    scale_fill_gradient2(limits = c(0, 1), low = "#000000", mid="red", high =  "#FFD700", na.value = "white", midpoint = 0.5) +
    geom_hline(
      aes(color = expected_motif, yintercept = 3.6),
      linewidth = 2
    ) +
    #scale_y_continuous(limits = c(0.5,1.6)) +
    facet_grid(~expected_motif, scales = "free_x", space = "free_x") +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.ticks.y = element_blank(),
      plot.title = element_markdown(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.background = element_blank(),   # Remove panel background,
      panel.margin = unit(0, "lines"),
      panel.border = element_rect(colour = "#141414", fill=NA, linewidth = 0),
      panel.spacing = unit(0, "lines"),  # Adjust spacing
      strip.placement = "bottom",
      strip.text.y.left = element_markdown(angle = 0, vjust = 0.5, hjust = 0, size =15),
      strip.background = element_blank(),
      strip.text.x.top = element_markdown()
    ) +
    labs(
      y = NULL,
      x = NULL,
      title = unique(df$name)
    ) +
    geom_text(data=df[identified == "No"], aes(label = "X", x=motif_markdown, y=1, color="X")) +
    geom_segment(data=df[identified == "No"], 
                 aes(x=motif_markdown, xend=motif_markdown, y=1.1, yend=1.4), 
                 color="white", size=2) +
    scale_color_manual(
      name = "", 
      breaks = c("Expected", "New",    "X"),
      values = c("Expected" = "#B0B0B0", "New" = "forestgreen",  "X" = "red"), 
      labels = c("Expected" = "Expected", "New" = "New", "X" = "Not identified"),
      limits = c("Expected", "New", "X"),
      drop = FALSE,
      guide = NULL
    )
  
  
  return(p)
}

shared_guide_color <- guides(
  color = guide_legend(
    override.aes = list(
      alpha = 1,
      shape = c(1, 1, 0),
      size = 5,
      color = c("#B0B0B0", "forestgreen", "red"),
      linetype = c(1, 1, 0),
      linewidth = c(1, 1, 0),
      label = c("", "", "x")
    )
  )
)

singles_plot <-  split(plot_data_separated, by = "name") %>% lapply(
  plot_single_isolate
)
p = wrap_plots(singles_plot, nrow = 2) + 
  shared_guide_color +
  plot_layout(guides = "collect")







gold_standard_motifs_scored <- gold_standard_motifs %>% 
  merge(
    bin_motifs_scored, 
    by.y = c(
      "bin",
      "motif",
      "mod_type",
      "mod_position"
    ),
    by.x = c(
      "sample_id",
      "motif",
      "mod_type",
      "mod_position"
    )
  )
gold_standard_motifs_scored[, group := NULL][, filepath := NULL][
  , mean := n_mod_bin / (n_mod_bin + n_nomod_bin)
]


# Modkit
monoculture_data_folder <- "data/monocultures"
monocultures_modkit <- list.files(monoculture_data_folder) %>% 
  lapply(function(x) {
    if (!file.exists(file.path(monoculture_data_folder, x))) return(NA)
    else {
      path = file.path(monoculture_data_folder, x, "modkit.motifs.tsv")
      if (!file.exists(path)) return(NA)
      motifs <- fread(path)
      motifs[
        , path := path
      ]
      
    }
  }) %>%  `[`(!is.na(.)) %>% rbindlist()
monocultures_modkit <- monocultures_modkit[
  , sample_type := "Isolates"
][
  , sample_id := str_extract(path, "monocultures/.+/") %>%  str_remove("monocultures") %>% str_remove_all("/")
][
  , identified := "Yes"
] %>% merge(
  gold_standard_motifs_scored,
  all = TRUE,
  by.x =c(
    "sample_id",
    "motif",
    "mod_code",
    "offset"
  ),
  by.y = c(
    "sample_id",
    "motif",
    "mod_type",
    "mod_position"
  )
)
monocultures_modkit[
  , identified := fifelse(
    is.na(identified),
    "No",
    "Yes"
  )
][
  is.na(gold_standard), gold_standard := FALSE
][
  , mean := fcase(
    is.na(mean), frac_mod,
    default = mean
  )
][
  , mod_type := mod_code
][
  , mod_type_pretty := map(mod_code, function(x) MOD_TYPE_PRETTY[[x]]) %>%  unlist()
][
  , motif_markdown := paste0(
    fifelse(offset <= 0, "", str_sub(motif, 1, offset)),
    "<strong>",
    str_sub(motif, offset+1, offset+1),
    "<sub><sup>",
    map(mod_code, function(x) MOD_TYPE_PRETTY[[x]]) %>%  unlist(),
    "</sup></sub>",
    "</strong>",
    str_sub(motif, offset+2, -1)
  )
][
  , name := str_remove(sample_id, "CVM\\d+_") %>% str_replace_all(
    "_", " "
  ) %>% str_remove(
    "DSMZ\\d+-"
  ) %>% str_replace(
    "e coli k12", "Escherichia coli MG1655"
  ) %>% str_replace( 
    "m ruber", "Meiothermus ruber"
  ) %>%  str_replace(
    "p thermoglucosidasius", "Parageobacillus thermoglucosidasius"
  ) %>% str_replace(
    "e coli dcm-dam", "Escherichia coli dam-/dcm-"
  ) %>% str_remove(
    " chromosome"
  ) %>% str_replace(
    "plasmidl", "plasmid1"
  ) %>% str_remove(
    " complete genome"
  ) %>% str_replace(
    "b2207", "B2207"
  ) %>% 
    str_remove(" Concatenation of \\d+ sequences") %>% 
    str_remove("final") %>% 
    str_replace_all("\\.+", " ") %>% 
    str_remove("^[A-Za-z]+ [A-Za-z]+-") %>% 
    str_remove("^[A-Za-z]+ [A-Za-z]+ [A-Z]+\\d+-")
][
  str_detect(name, "Escherichia coli$") & ((sample_type == "ZymoGut") | (sample_type == "ZymoHMW")), name := "Escherichia coli B1109"
][
  str_detect(name, "Escherichia coli plasmid$") & ((sample_type == "ZymoGut") | (sample_type == "ZymoHMW")), name := "Escherichia coli B1109 plasmid"
][
  str_detect(name, "Saccharomyces cerevisiae"), name := "Saccharomyces cerevisiae"
][
  str_detect(name, "coli"), name := str_replace(name, "E.coli", "Escherichia coli") 
][
  str_detect(name, "plasmid"), name := str_replace(name, "plasmid.*", "**\\0**")
][
  , pretty_modtype := fcase(
    mod_code == "a", "6mA",
    mod_code == "m", "5mC",
    mod_code == "21839", "4mC"
  )
][
  , name := str_replace(name, "[A-Za-z]+ [A-Za-z]+", "*\\0*")
][
  , gold_standard_chr := fifelse(gold_standard, "yes", "no") 
]
plot_data <- monocultures_modkit[name %in% isolate_selection][
  (gold_standard == TRUE) | (identified == "Yes")
][
  !str_detect(name, "natans") 
]

singles_plot <-  split(plot_data, by = "name") %>% lapply(
  plot_single_isolate
)
p = wrap_plots(singles_plot) + 
  shared_guide_color +
  plot_layout(guides = "collect")




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## Full ploy with all motifs
heatmap_plot <- function(df) {
  axis_orders <- get_axis_order(dcast(df, name ~ motif_markdown, value.var = "mean", fill = 0, fun.aggregate = mean))
  heatmap_plot <- ggplot(df) +
    aes(x = factor(motif_markdown, levels = axis_orders[["x"]]), y = factor(name, levels =  sort(axis_orders[["y"]], decreasing = TRUE))) +
    geom_tile(aes(fill = factor(identified, levels=c("Yes", "No"))), color = "#ffffff00", alpha = 0) +
    scale_fill_manual(values = c("No" = PLOT_COLORS[[4]], "Yes" = PLOT_COLORS[[2]])) +
    guides(fill = guide_legend(
      title = "Identified",
      override.aes = list(alpha = 1))
    ) +
    new_scale("fill") +
    
    geom_tile(aes(fill = mean), linewidth = 0, alpha = 0, color = "#33333300") +
    scale_fill_gradient2(
      limits = c(0, 1), low = "white", high = "gray20",
      labels = scales::percent
    ) +
    new_scale("fill") +
    
    geom_tile(aes(fill = mean),  color = "black",linewidth = 0.1) +
    scale_fill_gradient2(guide = "none", limits = c(0, 1), low = "white", high =  PLOT_COLORS[4], na.value = "white") +
    new_scale("fill") +
    
    geom_tile(data=df[identified == "Yes"], aes(fill = mean),  color = "black",linewidth = 0.1) + 
    scale_fill_gradient2(guide = "none", limits = c(0, 1), low = "white", high =  PLOT_COLORS[2], na.value = "white") +
    
    geom_text(data=df[(n_mod_bin.scored + n_nomod_bin.scored)==0], label = ".", hjust=0.5, vjust=0, size=0.5) +
    facet_grid(sample_type~pretty_modtype, space="free", scale="free") +
    theme(
      axis.text.y = element_markdown(size = 6),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.background = element_blank(),   # Remove panel background,
      panel.margin = unit(0.3, "lines"),
      panel.border = element_rect(colour = "#141414", fill=NA, linewidth = 0),
      strip.placement = "outside",
      strip.text.y.left = element_markdown(angle = 0, vjust = 0.5, hjust = 0, size =15),
      strip.background = element_blank()
    ) +
    labs(
      y = NULL,
      x = "Motifs"
    ) +
    geom_text(aes(color = "No motif observations"), alpha = 0, label ="") +
    geom_text(aes(color = "Assigned MTase"), alpha=0, label = "") +
    geom_text(data=df[(n_mod_bin.scored + n_nomod_bin.scored)==0], label = ".", hjust=0.5, vjust=0.5, size=3) +
    geom_text(data=df[linked == TRUE], label = "+", col = "orange", hjust=0.5, vjust=0.5, size=3) +
    scale_color_manual(name = "Indicator", values = c("No motif observations" = "black", "Assigned MTase" = "orange")) +
    guides(color = guide_legend(override.aes = list(alpha = 1, label = c("+", "●"), size = c(3, 2),  fill = c("orange", "black")))) 
  return(heatmap_plot)
}
heatmap_plot_goldstandard <- function(df) {
  axis_orders <- get_axis_order(dcast(df, name ~ motif_markdown, value.var = "mean", fill = 0, fun.aggregate = mean))
  heatmap_plot <- ggplot(df) +
    aes(x = factor(motif_markdown, levels = axis_orders[["x"]]), y = factor(name, levels =  sort(axis_orders[["y"]], decreasing = TRUE))) +
    geom_tile(aes(fill = factor(identified, levels=c("Yes", "No"))), color = "#ffffff00", alpha = 0) +
    scale_fill_manual(values = c("No" = PLOT_COLORS[[4]], "Yes" = PLOT_COLORS[[2]])) +
    guides(fill = guide_legend(
      title = "Identified",
      override.aes = list(alpha = 1))
    ) +
    new_scale("fill") +
    
    geom_tile(aes(fill = mean, color = gold_standard), linewidth = 0, alpha = 0, color = "#33333300") +
    scale_fill_gradient2(
      limits = c(0, 1), low = "white", high = "gray20",
      labels = scales::percent
    ) +
    new_scale("fill") +
    geom_tile(aes(fill = mean, color = gold_standard), linewidth = 0.5) +
    scale_fill_gradient2(guide = "none", limits = c(0, 1), low = "white", high =  PLOT_COLORS[4], na.value = "white") +
    new_scale("fill") +
    geom_tile(data=df[identified == "Yes"], aes(fill = mean, color = gold_standard),  linewidth = 0.5) + 
    scale_fill_gradient2(guide = "none", limits = c(0, 1), low = "white", high =  PLOT_COLORS[2], na.value = "white") +
    scale_color_manual(name = "Gold standard", values = c("gray60", "orange")) +
    new_scale("color") +
    geom_text(data=df[(n_mod_bin.scored + n_nomod_bin.scored)==0], label = ".", hjust=0.5, vjust=0, size=0.5) +
    facet_grid(sample_type~pretty_modtype, space="free", scale="free") +
    theme(
      axis.text.y = element_markdown(size = 9),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.background = element_blank(),   # Remove panel background,
      panel.margin = unit(0.3, "lines"),
      panel.border = element_rect(colour = "#141414", fill=NA, linewidth = 0),
      strip.placement = "outside",
      strip.text.y.left = element_markdown(angle = 0, vjust = 0.5, hjust = 0, size =15),
      strip.background = element_blank()
    ) +
    labs(
      y = NULL,
      x = "Motifs"
    ) 
  return(heatmap_plot)
}



plot_data <- bin_motifs_full
heatmap_full <- heatmap_plot(plot_data)

ggsave(
  "figures/motif_heatmap.png",
  heatmap_full,
  width = 21,
  height = 8
)


################################################################################
# Subset e coli
plot_data_ecoli <- plot_data[TRUE][
  str_detect(bin, "coli")
][
  , keep := any(identified == "Yes"), by = .(motif_markdown)
][keep == TRUE][, keep := NULL]
heatmap_ecoli <- heatmap_plot(plot_data_ecoli)
ggsave(
  "figures/motif_heatmap_ecoli.png",
  heatmap_ecoli,
  width = 9,
  height = 5
)

# subset gold standard
gold_standard <- c("Zymomonas mobilis", "Cellulophaga lytica", "Escherichia coli", "Kangiella aquimarina", "Meiothermus ruber", "Shewanella oneidensis")
plot_data_gold_standard <- plot_data[TRUE][
  sample_type == "Isolates"
][
  , keep := any(gold_standard == TRUE), by = .(name)
][keep == TRUE][, keep := NULL][
  , keep := any(identified == "Yes") | any(gold_standard == TRUE), by = .(motif_markdown)
][keep == TRUE][, keep := NULL]
heatmap_gold_standard <- heatmap_plot_goldstandard(plot_data_gold_standard)
ggsave(
  "figures/motif_heatmap_gold_standard.png",
  heatmap_gold_standard,
  width = 17,
  height = 5
)

# Isolates
plot_data_isolates <- plot_data[TRUE][
  sample_type == "Isolates"
][
  , keep := any(identified == "Yes"), by = .(motif_markdown)
][keep == TRUE][, keep := NULL]
heatmap_isolates <- heatmap_plot(plot_data_isolates)
ggsave(
  "figures/motif_heatmap_isolates.png",
  heatmap_isolates,
  width = 15,
  height = 6
)



