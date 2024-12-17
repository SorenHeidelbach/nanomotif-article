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
  "rjson",
  "gt",
  "gtExtras"
)
getwd()
source("src/constants.R")
source("src/themes.R")
source("src/utility.R")

# Load contig bin from sample name
load_contig_bin <- function(sample) {
  path <- paste0("data/real_communities/", sample, "/mmlong2_lite/tmp/binning/contig_bin.tsv")
  contig_bin <- fread(path, header = FALSE)
  setnames(contig_bin, c("contig", "bin"))
  return(contig_bin)
}

load_bin_info <- function(path) {
  bin_info <- fread(path)
  # Ensure columns are present:
  # If contamination_checkm2 & completeness_checkm2 exist, use them.
  # Else assume Completeness & Contamination exist.
  if("contamination_checkm2" %in% colnames(bin_info)){
    bin_info[
      , bin_quality := case_when(
        completeness_checkm2 > 90 & contamination_checkm2 < 5 ~ "HQ",
        completeness_checkm2 > 50 & contamination_checkm2 < 10 ~ "MQ",
        TRUE ~ "LQ"
      )
    ]
  } else {
    bin_info[
      , bin_quality := case_when(
        Completeness > 90 & Contamination < 5 ~ "HQ",
        Completeness > 50 & Contamination < 10 ~ "MQ",
        TRUE ~ "LQ"
      )
    ]
  }
  return(bin_info)
}

load_bin_motifs_w_info <- function(
    motifs_bin_path,
    bin_info_path
){
  bin_info <- load_bin_info(bin_info_path)
  bin_motifs <- fread(motifs_bin_path) 
  output <- bin_motifs %>% merge(
    bin_info,
    by = "bin",
    all.y = TRUE
  )
  return(output)
}

load_sample <- function(sample){
  bin_tsv_path <- fifelse(
    sample == "fecal_simple_recovered", "mmlong2_bins.tsv", "mmlong2_lite_bins.tsv"
  )
  version <- fifelse(
    sample == "mfd02199_backup", "0.4.13", "0.4.16"
  )
  sample_data <- load_bin_motifs_w_info(
    paste0("../data/real_communities/", sample, "/nanomotif_", version,"/bin-motifs.tsv"),
    paste0("../data/real_communities/", sample, "/mmlong2_lite/results/", bin_tsv_path)
  )
  return(sample_data)
}

load_tiara_filt_contigs <- function(tiara_path) {
  # Ensure tiara_path file has one column
  df <- fread(tiara_path, header = FALSE)
  names(df) <- c("contig")
  return(df$contig)
}

load_valid_contigs <- function(mmlong2_cov_all_path, mmlong2_tiara_path, min_length=3000, min_cov=0) {
  tiara_filt_contigs <- load_tiara_filt_contigs(mmlong2_tiara_path)
  # Ensure cov_all.tsv has columns: contig, length, cov, cov_np, cov_var
  df <- fread(mmlong2_cov_all_path)
  names(df) <- c("contig", "length", "cov", "cov_np", "cov_var")
  df <- df[
    length > min_length
  ][
    cov > min_cov
  ][
    contig %in% tiara_filt_contigs
  ]
  return(df)
}

load_read_stats <- function(path, sample) {
  # Ensure file at path has two columns: key and value.
  read_stats <- fread(path, header = FALSE)[, Sample := sample] %>% dcast(Sample ~ V1, value.var = "V2")
  return(read_stats)
}

load_binnary_contaminants <- function(
    binnary_contaminants_path,
    bin_info_path,
    sample,
    min_cov = 10
) {
  # Ensure this file has a 'contig' and 'bin' column.
  df <- fread(binnary_contaminants_path)
  
  contaminants <- unique(df$contig)
  bin_info <- load_bin_info(bin_info_path)
  
  hq_bins <- bin_info[bin_quality == "HQ" & (cov >= min_cov)]$bin
  hq_contigs <- df[bin %in% hq_bins]$contig
  contaminants_hq <- contaminants[contaminants %in% hq_contigs]
  
  mq_bins <- bin_info[bin_quality == "MQ" & (cov >= min_cov)]$bin
  mq_contigs <- df[bin %in% mq_bins]$contig
  contaminants_mq <- contaminants[contaminants %in% mq_contigs]
  
  contaminants_low_cov <- contaminants[!(contaminants %in% mq_contigs) & !(contaminants %in% hq_contigs)]
  
  median_df <- df[!duplicated(contig)][
    , .(
      n_contaminant = .N
    ), by = bin
  ]
  data.table(
    sample = sample,
    contaminant_count = length(unique(contaminants)),
    contaminant_count_in_hq = length(contaminants_hq),
    contaminant_count_in_mq = length(contaminants_mq),
    contaminant_count_in_low_cov = length(contaminants_low_cov),
    median_contaminant_hq = median(median_df[bin %in% hq_bins]$n_contaminant),
    median_contaminant_mq = median(median_df[bin %in% mq_bins]$n_contaminant)
  )
}

load_genomad_classification <- function(
    genomad_path,
    threshold
) {
  # Ensure genomad file has columns: plasmid_score, virus_score, chromosome_score, seq_name
  genomad <- fread(genomad_path)[
    , classification := fcase(
      plasmid_score > threshold, "plasmid",
      virus_score > threshold, "virus",
      chromosome_score > threshold, "chromosome",
      default = "none"
    )
  ]
  return(genomad)
}

load_unbinned_plasmid_contigs <- function(
    tiara_path,
    mmlong2_cov_all_path,
    genomad_path,
    contig_bin_path,
    threshold, 
    min_cov = 10, 
    min_length = 10000
) {
  genomad <- load_genomad_classification(genomad_path, threshold)
  valid_contigs <- load_valid_contigs(
    mmlong2_tiara_path = tiara_path,
    mmlong2_cov_all_path = mmlong2_cov_all_path,
    min_cov = min_cov,
    min_length = min_length
  )
  
  contig_bin <- fread(contig_bin_path, header = FALSE)
  # Contig_bin assumed to have two columns: contig and bin
  plasmid_contigs <- genomad[classification %in% c("plasmid", "virus")]$seq_name
  # Filter out contigs that are already binned
  plasmid_contigs <- plasmid_contigs[! (plasmid_contigs %in% contig_bin[[1]])]
  # Filter only those that appear in valid_contigs$contig
  plasmid_contigs <- plasmid_contigs[plasmid_contigs %in% valid_contigs$contig]
  return(plasmid_contigs)
}

load_binnary_include <- function(
    binnary_include_path, 
    tiara_path,
    mmlong2_cov_all_path,
    bin_info_path,
    genomad_path,
    contig_bin_path,
    sample,
    min_cov=10
) {
  # Ensure include_contigs.tsv has columns: contig, assigned_bin, confidence
  included <- fread(binnary_include_path)[confidence == "high_confidence"]
  bin_info <- load_bin_info(bin_info_path)
  
  hq_bins <- bin_info[bin_quality == "HQ" & (cov >= min_cov)]$bin
  included_in_hq <- included[assigned_bin %in% hq_bins]
  
  mq_bins <- bin_info[bin_quality == "MQ" & (cov >= min_cov)]$bin
  included_in_mq <- included[assigned_bin %in% mq_bins]
  
  unbinned_plasmid_contigs <- load_unbinned_plasmid_contigs(
    tiara_path,
    mmlong2_cov_all_path,
    genomad_path,
    contig_bin_path,
    threshold = 0.75
  )
  included_plasmids <- included$contig[included$contig %in% unbinned_plasmid_contigs]
  
  median_df <- included[!duplicated(contig)][
    , .(
      n_included = .N,
      n_included_mge = length(unique(contig[contig %in% included_plasmids]))
    ), by = assigned_bin
  ]
  data.table(
    sample = sample,
    inclusion_count = length(unique(included$contig)),
    unbinned_plasmid_contigs = length(unbinned_plasmid_contigs),
    unbinned_plasmid_contigs_included = length(unique(included_plasmids)),
    n_included_in_hq = length(unique(included_in_hq$contig)),
    n_included_in_mq = length(unique(included_in_mq$contig)),
    median_included_hq = median(median_df[assigned_bin %in% hq_bins]$n_included),
    median_included_mq = median(median_df[assigned_bin %in% mq_bins]$n_included),
    median_included_mge = median(median_df[(assigned_bin %in% hq_bins) | assigned_bin %in% mq_bins]$n_included_mge),
    median_included_mge_hq = median(median_df[assigned_bin %in% hq_bins]$n_included_mge),
    median_included_mge_mq = median(median_df[assigned_bin %in% mq_bins]$n_included_mge)
    
  )
}

get_mapped_bases <- function(bin_info) {
  if ("genome_size" %in% colnames(bin_info)) {
    bin_info[
      , .(mapped_bases = sum(cov * genome_size))
    ]$mapped_bases
  } else {
    bin_info[
      , .(mapped_bases = sum(cov * Genome_Size))
    ]$mapped_bases
  }
}

get_assembly_stats <- function(
    mmlong2_cov_all_path, 
    mmlong2_tiara_path, 
    binnary_contaminants_path,
    binnary_include_path,
    genomad_path,
    bin_info_path,
    contig_bin_path,
    sample
){
  valid_contigs <- load_valid_contigs(mmlong2_cov_all_path, mmlong2_tiara_path)
  contig_bin <- load_contig_bin(sample)
  binned_contigs <- contig_bin$contig
  bin_info <- load_bin_info(bin_info_path)
  
  binned_bases <- get_mapped_bases(bin_info)
  # General assembly stats
  assembly_stats <- data.table(
    n_contigs = length(valid_contigs$contig),
    n_contig_binned = sum(valid_contigs$contig %in% binned_contigs),
    n50_contigs = N50(valid_contigs$length),
    sample = sample
  )
  
  # Contaminant
  contamination <- load_binnary_contaminants(
    binnary_contaminants_path,
    bin_info_path,
    sample = sample
  )
  
  # Inclusion 
  inclusion <- load_binnary_include(
    binnary_include_path, 
    mmlong2_tiara_path,
    mmlong2_cov_all_path,
    bin_info_path,
    genomad_path,
    contig_bin_path,
    sample = sample
  )
  collected_output <- merge(
    contamination,
    inclusion,
    by = "sample"
  ) %>%  merge(
    assembly_stats, by = "sample")
  collected_output[
    , binned_bases := binned_bases
  ]
  return(collected_output)
}




samples <- c(
  "Fecal, Simple"="fecal_simple_recovered",
  "Fecal, Complex"="fecal_inhouse",
  "ZymoFecal" = "ZymoFecal",
  "Anaerobic digester"="anaerobic_digester",
  "Soil"="mfd02199_backup"
)


info <- rbind(
  get_assembly_stats(
    mmlong2_cov_all_path = "../data/real_communities/fecal_simple_recovered/mmlong2_lite/tmp/binning/mapping/cov_all.tsv", 
    mmlong2_tiara_path = "../data/real_communities/fecal_simple_recovered/mmlong2_lite/tmp/filtering/contigs_filt_len_euk.txt", 
    binnary_contaminants_path = "../data/real_communities/fecal_simple_recovered/nanomotif_binnary/bin_contamination.tsv",
    binnary_include_path = "../data/real_communities/fecal_simple_recovered/nanomotif_binnary/include_contigs.tsv",
    genomad_path = "../data/real_communities/fecal_simple_recovered/genomad/mmlong2_assembly_aggregated_classification/mmlong2_assembly_aggregated_classification.tsv",
    bin_info_path = "../data/real_communities/fecal_simple_recovered/mmlong2_lite/results/mmlong2_bins.tsv",
    contig_bin_path = "data/real_communities/fecal_simple_recovered/mmlong2_lite/tmp/binning/contig_bin.tsv",
    sample = "fecal_simple_recovered"
  ),
  get_assembly_stats(
    mmlong2_cov_all_path = "../data/real_communities/ZymoFecal/mmlong2_lite/tmp/binning/mapping/cov_all.tsv", 
    mmlong2_tiara_path = "../data/real_communities/ZymoFecal/mmlong2_lite/tmp/filtering/contigs_filt_len_euk.txt", 
    binnary_contaminants_path = "../data/real_communities/ZymoFecal/nanomotif_binnary/bin_contamination.tsv",
    binnary_include_path = "../data/real_communities/ZymoFecal/nanomotif_binnary/include_contigs.tsv",
    genomad_path = "../data/real_communities/ZymoFecal/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv",
    bin_info_path = "../data/real_communities/ZymoFecal/mmlong2_lite/results/mmlong2_lite_bins.tsv",
    contig_bin_path = "data/real_communities/ZymoFecal/mmlong2_lite/tmp/binning/contig_bin.tsv",
    sample = "ZymoFecal"
  ),
  get_assembly_stats(
    mmlong2_cov_all_path = "../data/real_communities/fecal_inhouse/mmlong2_lite/tmp/binning/mapping/cov_all.tsv", 
    mmlong2_tiara_path = "../data/real_communities/fecal_inhouse/mmlong2_lite/tmp/filtering/contigs_filt_len_euk.txt", 
    binnary_contaminants_path = "../data/real_communities/fecal_inhouse/nanomotif_binnary/bin_contamination.tsv",
    binnary_include_path = "../data/real_communities/fecal_inhouse/nanomotif_binnary/include_contigs.tsv",
    genomad_path = "../data/real_communities/fecal_inhouse/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv",
    bin_info_path = "../data/real_communities/fecal_inhouse/mmlong2_lite/results/mmlong2_lite_bins.tsv",
    contig_bin_path = "data/real_communities/fecal_inhouse/mmlong2_lite/tmp/binning/contig_bin.tsv",
    sample = "fecal_inhouse"
  ),
  get_assembly_stats(
    mmlong2_cov_all_path = "../data/real_communities/anaerobic_digester/mmlong2_lite/tmp/binning/mapping/cov_all.tsv", 
    mmlong2_tiara_path = "../data/real_communities/anaerobic_digester/mmlong2_lite/tmp/filtering/contigs_filt_len_euk.txt", 
    binnary_contaminants_path = "../data/real_communities/anaerobic_digester/nanomotif_binnary/bin_contamination.tsv",
    binnary_include_path = "../data/real_communities/anaerobic_digester/nanomotif_binnary/include_contigs.tsv",
    genomad_path = "../data/real_communities/anaerobic_digester/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv",
    bin_info_path = "../data/real_communities/anaerobic_digester/mmlong2_lite/results/mmlong2_lite_bins.tsv",
    contig_bin_path = "data/real_communities/anaerobic_digester/mmlong2_lite/tmp/binning/contig_bin.tsv",
    sample = "anaerobic_digester"
  ),
  get_assembly_stats(
    mmlong2_cov_all_path = "../data/real_communities/mfd02199_backup/mmlong2_lite/tmp/binning/cov_tmp/1_cov.tsv", 
    mmlong2_tiara_path = "../data/real_communities/mfd02199_backup/mmlong2_lite/tmp/eukfilt/contigs_filt.txt", 
    binnary_contaminants_path = "../data/real_communities/mfd02199_backup/nanomotif_binnary/bin_contamination.tsv",
    binnary_include_path = "../data/real_communities/mfd02199_backup/nanomotif_binnary/include_contigs.tsv",
    genomad_path = "../data/real_communities/mfd02199_backup/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv",
    bin_info_path = "../data/real_communities/mfd02199_backup/mmlong2_lite/results/mmlong2_lite_bins.tsv",
    contig_bin_path = "data/real_communities/mfd02199_backup/mmlong2_lite/tmp/binning/contig_bin.tsv",
    sample = "mfd02199_backup"
  )
)

read_info <- rbind(
  load_read_stats("data/real_communities/fecal_simple/qc_reads/nanostat", "fecal_simple"),
  load_read_stats("data/real_communities/ZymoFecal/qc_reads/nanostat", "ZymoFecal"),
  load_read_stats("data/real_communities/fecal_inhouse/qc_reads/nanostat", "fecal_inhouse"),
  load_read_stats("data/real_communities/anaerobic_digester/qc_reads/nanostat", "anaerobic_digester")
)

motif_info_bin_level <- names(samples) %>% lapply(function(x) {
  df <- load_sample(samples[x])
  df[, sample := samples[x]]
  df <- df[
    cov >= 10
  ]
  df_median <- df[!is.na(motif)][
    , .(
      n_motif = length(motif[!is.na(motif)]),
      n_motif_hq = length(motif[!is.na(motif) & (bin_quality == "HQ")]),
      n_motif_mq = length(motif[!is.na(motif) & (bin_quality == "MQ")]),
      n_motif_5mC = length(motif[!is.na(motif) & (mod_type == "m")]),
      n_motif_4mC = length(motif[!is.na(motif) & (mod_type == "21839")]),
      n_motif_6mA = length(motif[!is.na(motif) & (mod_type == "a")])
    ), by = bin
  ]
  return(df_median)
}) %>% rbindlist()
motif_info_bin_level[][
  , .(
    median_value = median(n_motif),
    avg = mean(n_motif)
  )
]
motif_info <- names(samples) %>% lapply(function(x) {
  df <- load_sample(samples[x])
  df[, sample := samples[x]]
  df <- df[
    cov >= 10
  ]
  df_median <- df[!is.na(motif)][
    , .(
      n_motif = length(motif[!is.na(motif)]),
      n_motif_hq = length(motif[!is.na(motif) & (bin_quality == "HQ")]),
      n_motif_mq = length(motif[!is.na(motif) & (bin_quality == "MQ")]),
      n_motif_5mC = length(motif[!is.na(motif) & (mod_type == "m")]),
      n_motif_4mC = length(motif[!is.na(motif) & (mod_type == "21839")]),
      n_motif_6mA = length(motif[!is.na(motif) & (mod_type == "a")])
    ), by = bin
  ][
    , .(
      n_motif_median =     median(n_motif[n_motif != 0]),
      n_motif_hq_median =  median(n_motif_hq[n_motif_hq != 0]),
      n_motif_mq_median =  median(n_motif_mq[n_motif_mq != 0]),
      n_motif_5mC_median = median(n_motif_5mC[n_motif != 0]),
      n_motif_4mC_median = median(n_motif_4mC[n_motif != 0]),
      n_motif_6mA_median = median(n_motif_6mA[n_motif != 0])
    )
  ][, sample := samples[x]]
  df_sum <- df[
    , .(
      sample = unique(sample),
      n_bin = unique(bin) %>%  length(),
      n_bin_hq = unique(bin[bin_quality == "HQ"]) %>%  length(),
      n_bin_mq = unique(bin[bin_quality == "MQ"]) %>%  length(),
      n_bin_hq_with_motif = unique(bin[!is.na(motif) & (bin_quality == "HQ")]) %>%  length(),
      n_bin_mq_with_motif = unique(bin[!is.na(motif) & (bin_quality == "MQ")]) %>%  length(),
      n_motif = length(motif[!is.na(motif)]),
      n_motif_hq = length(motif[!is.na(motif) & (bin_quality == "HQ")]),
      n_motif_mq = length(motif[!is.na(motif) & (bin_quality == "MQ")]),
      n_motif_5mC = length(motif[!is.na(motif) & (mod_type == "m")]),
      n_motif_4mC = length(motif[!is.na(motif) & (mod_type == "21839")]),
      n_motif_6mA = length(motif[!is.na(motif) & (mod_type == "a")])
      
    )
  ][
    , n_motif_pr_hq_bin := n_motif_hq / n_bin_hq
  ][
    , n_motif_pr_mq_bin := n_motif_mq / n_bin_mq
  ]
  df_out <- merge(
    df_sum,
    df_median
  )
  return(df_out)
}) %>%  rbindlist()


full_summary <- merge(
  motif_info, 
  read_info,
  by.x = "sample",
  by.y = "Sample",
  all.x = TRUE
) %>% merge(
  info,
  by.x = "sample",
  by.y = "sample"
)















full_summary[
  , n_bins_formatted := paste0(n_bin_mq_with_motif, " | ", n_bin_hq_with_motif)
][
  , n_bin_with_motif_formatted := paste0(n_bin_mq_with_motif, " | ", n_bin_hq_with_motif)
][
  , percent_mag_with_motif_formatted := paste0(round(100*n_bin_mq_with_motif/n_bin_mq, 0), "% | ", round(100*n_bin_hq_with_motif/n_bin_hq, 0), "%")
][
  , n_motif_hq_included := paste0(n_motif, " (", n_motif_hq, ")")
][
  , gb_binned := binned_bases / 1e9
][
  , gb_sequenced := paste0(round(as.numeric(number_of_bases)/1e9, 1), " (", round(100*binned_bases/as.numeric(number_of_bases), 0), "%)")
][
  , n50_kbp := round(n50_contigs/1000, 0)
][
  , motif_pr_bin_formatted := paste0(round(n_motif_pr_mq_bin, 2), " | ", round(n_motif_pr_hq_bin, 2))
][
  , n_5mC_pr_bin := round(n_motif_5mC / n_bin, 2)
][
  , n_4mC_pr_bin := round(n_motif_4mC / n_bin, 2)
][
  , n_6mA_pr_bin := round(n_motif_6mA / n_bin, 2)
][
  , n_type_pr_bin := paste(n_6mA_pr_bin, n_5mC_pr_bin, n_4mC_pr_bin, sep = " | ")
][
  , contaminant_formatted := paste0(contaminant_count_in_mq, " | ", contaminant_count_in_hq)
][
  , inclusion_formatted := paste0(n_included_in_mq, " | ", n_included_in_hq)
][
  , mge_formatted := paste0(unbinned_plasmid_contigs_included, " (", round(100*unbinned_plasmid_contigs_included/unbinned_plasmid_contigs, 0), "%)")
][
  , contigs_formatted := paste0(n_contigs, " (", round(100*n_contig_binned/n_contigs, 0),  "%)")
]

data_to_present <- c(
  "Sample"="sample",
  "binned_bases"="gb_binned",
  "# Bases [Gbp]"="gb_sequenced",
  "# Contigs (%binned)"="contigs_formatted",
  "N50 [kbp]"="n50_kbp",
  "# Bins MQ|HQ"="n_bins_formatted",
  "% Bin w\\ motif MQ|HQ"="percent_mag_with_motif_formatted",
  "Motifs/bin MQ|HQ"="motif_pr_bin_formatted",
  "Motif pr. bin 6mA|5mC|4mC"="n_type_pr_bin",
  "# Contaminants removed MQ|HQ"="contaminant_formatted",
  "# Contigs included MQ|HQ"="inclusion_formatted",
  "MGEs included (% of all MGEs)"="mge_formatted"
)
summary_table <- full_summary %>% select(data_to_present) %>%  gt() %>%  
  cols_align('right') %>% 
  cols_label(.fn = html)
gtsave(
  summary_table,
  "figures/metagenome_summary_table.html"
)
 
















plasmid_contigs <- rbind(
  data.table(
    sample = "fecal_simple_recovered",
    mge_contigs = load_unbinned_plasmid_contigs(
      "../data/real_communities/fecal_simple_recovered/mmlong2_lite/tmp/filtering/contigs_filt_len_euk.txt", 
      mmlong2_cov_all_path = "../data/real_communities/fecal_simple_recovered/mmlong2_lite/tmp/binning/mapping/cov_all.tsv", 
      genomad_path = "../data/real_communities/fecal_simple_recovered/genomad/mmlong2_assembly_aggregated_classification/mmlong2_assembly_aggregated_classification.tsv",
      contig_bin_path = "data/real_communities/fecal_simple_recovered/mmlong2_lite/tmp/binning/contig_bin.tsv",
      0.75,
      min_cov = 10,
      min_length = 10000
      )),
    data.table(
      sample = "ZymoFecal",
      mge_contigs = load_unbinned_plasmid_contigs(
      "../data/real_communities/ZymoFecal/mmlong2_lite/tmp/filtering/contigs_filt_len_euk.txt", 
      mmlong2_cov_all_path = "../data/real_communities/ZymoFecal/mmlong2_lite/tmp/binning/mapping/cov_all.tsv", 
      genomad_path = "../data/real_communities/ZymoFecal/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv",
      contig_bin_path = "data/real_communities/ZymoFecal/mmlong2_lite/tmp/binning/contig_bin.tsv",
      0.75,
      min_cov = 10,
      min_length = 10000
      )),
    data.table(
      sample = "fecal_inhouse",
      mge_contigs = load_unbinned_plasmid_contigs(
      "../data/real_communities/fecal_inhouse/mmlong2_lite/tmp/filtering/contigs_filt_len_euk.txt", 
      mmlong2_cov_all_path = "../data/real_communities/fecal_inhouse/mmlong2_lite/tmp/binning/mapping/cov_all.tsv", 
      genomad_path = "../data/real_communities/fecal_inhouse/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv",
      contig_bin_path = "data/real_communities/fecal_inhouse/mmlong2_lite/tmp/binning/contig_bin.tsv",
      0.75,
      min_cov = 10,
      min_length = 10000
      )),
    data.table(
      sample = "anaerobic_digester",
      mge_contigs = load_unbinned_plasmid_contigs(
      "../data/real_communities/anaerobic_digester/mmlong2_lite/tmp/filtering/contigs_filt_len_euk.txt", 
      mmlong2_cov_all_path = "../data/real_communities/anaerobic_digester/mmlong2_lite/tmp/binning/mapping/cov_all.tsv", 
      genomad_path = "../data/real_communities/anaerobic_digester/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv",
      contig_bin_path = "data/real_communities/anaerobic_digester/mmlong2_lite/tmp/binning/contig_bin.tsv",
      0.75,
      min_cov = 10,
      min_length = 10000
      )),
    data.table(
      sample = "mfd02199_backup",
      mge_contigs = load_unbinned_plasmid_contigs(
      "../data/real_communities/mfd02199_backup/mmlong2_lite/tmp/eukfilt/contigs_filt.txt", 
      mmlong2_cov_all_path = "../data/real_communities/mfd02199_backup/mmlong2_lite/tmp/binning/cov_tmp/1_cov.tsv", 
      genomad_path = "../data/real_communities/mfd02199_backup/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv",
      contig_bin_path = "data/real_communities/mfd02199_backup/mmlong2_lite/tmp/binning/contig_bin.tsv",
      0.75,
      min_cov = 10,
      min_length = 10000
      ))
  )




bin_sup_df <- names(samples) %>% lapply(function(x) {
  df <- load_sample(samples[x])
  names(df)[base::grepl("enome_.ize", names(df))] <- c("genome_size")
  if (any(base::grepl("GC_Content", names(df)))) {
    df[
      , GC_Content := GC_Content * 100
    ]
  }
  names(df)[base::grepl("[Gg][Cc]", names(df))] <- c("gc")
  names(df)[base::grepl(".ontig_.50", names(df))] <- c("contig_n50")
  names(df)[base::grepl("Completeness", names(df)) | base::grepl("completeness_checkm2", names(df))] <- c("completeness")
  names(df)[base::grepl("Contamination", names(df)) | base::grepl("contamination_checkm2", names(df))] <- c("contamination")
  
  df <- df[!is.na(motif)][, sample := x][
      , .(
        n_motif = length(motif[!is.na(motif)]),
        n_motif_5mC = length(motif[!is.na(motif) & (mod_type == "m")]),
        n_motif_4mC = length(motif[!is.na(motif) & (mod_type == "21839")]),
        n_motif_6mA = length(motif[!is.na(motif) & (mod_type == "a")]),
        cov= unique(cov),
        genome_size = unique(genome_size),
        gc = unique(gc),
        contig_n50 = unique(contig_n50),
        completeness = unique(completeness),
        contamination = unique(contamination)
      ), by = .(bin, sample)
  ]

  path <- paste0("data/real_communities/", samples[x], "/nanomotif_binnary/include_contigs.tsv")
  included <- fread(path)[confidence == "high_confidence"]
  include <- included[!duplicated(contig)][
    , .(
      n_contig_included = .N,
      n_mge_contigs_included = sum(contig %in% plasmid_contigs[sample == samples[x]]$mge_contigs)
    ), by = assigned_bin
  ][, sample := samples[x]][, bin := assigned_bin][, assigned_bin := NULL]
  
  df <- merge(
    df, 
    include, all.x = TRUE
  )[
    is.na(n_contig_included), n_contig_included := 0
  ][
    is.na(n_mge_contigs_included), n_mge_contigs_included := 0
  ]
  contamination_path <- paste0("data/real_communities/", samples[x], "/nanomotif_binnary/bin_contamination.tsv")
  contamination <- fread(contamination_path)[!duplicated(contig)][, .(n_contigs_removed = .N), by= bin]
  print(contamination)
  df <- merge(
    df, 
    contamination, all.x = TRUE, by = "bin"
  )[
    is.na(n_contigs_removed), n_contigs_removed := 0
  ]
  
  return(df)
}) %>%  rbindlist()



fwrite(
  bin_sup_df,
  file = "supplementary_data_2.tsv", sep = "\t"
)
fwrite(
  plasmid_contigs[
    , sample := unique(names(samples[samples == sample])), by = sample
  ],
  "supplementary_data_3_mges.tsv", sep = "\t"
)
