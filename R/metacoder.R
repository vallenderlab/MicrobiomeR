


get_metacoder_obj <- function(phyloseq_object=NULL, rdata_file=NULL, filter_flag = FALSE,
                              ambiguous_filter=c("uncultured", "Unassigned", "Ambiguous"),
                              kingdom_filter=c("Bacteria")) {
  if (!is.null(rdata_file)){
    # Load data
    load(rdata_file)
    # Reassign loaded metacoder object to "master_metacoder"
    metacoder_obj <- get("metacoder_obj")
  } else {
    # Convert Phyloseq objects to Metacoder objects
    metacoder_obj <- metacoder::parse_phyloseq(phyloseq_object)
    # Calculate a proportionate abundance value for each sample
    metacoder_obj$data$otu_table <- metacoder::calc_obs_props(metacoder_obj, dataset = "otu_table", cols = metacoder_obj$data$sam_data$sample_ids, other_cols = TRUE)
    # Calculate the proportionate abundace per-taxon value for each sample
    metacoder_obj$data$tax_table <- metacoder::calc_taxon_abund(metacoder_obj, dataset = "otu_table", cols = metacoder_obj$data$sam_data$sample_ids)
    # Calculate the differences between Stress vs Control
    metacoder_obj$data$diff_table <- metacoder::compare_groups(metacoder_obj, dataset = "tax_table", cols = metacoder_obj$data$sam_data$sample_ids, groups = metacoder_obj$data$sam_data$TreatmentGroup, func = new_comparison, other_cols = TRUE)
    # Rename the tax_table columns for future filtering or manipulations
    colnames(metacoder_obj$data$tax_table) <- sapply(colnames(metacoder_obj$data$tax_table), function(x) if(x!="taxon_id"){sprintf("T%s", x)}else{x})
  }

  if (filter_flag == TRUE) {
    # Put data through initial filtering of ambiguous and kingdom data
    # Note always do this after calculating relative/proportionate abundance
    processed_metacoder <- initial_filter(
      obj = metacoder_obj,
      ambiguous_filter = ambiguous_filter,
      kingdom_filter = kingdom_filter
    )
    return(processed_metacoder)
  } else {
    return(metacoder_obj)
  }
}



new_comparison <- function(abund_1, abund_2) {
  log_med_ratio <- log2(median(abund_1) / median(abund_2))
  if (is.nan(log_med_ratio)) {
    log_med_ratio <- 0
  }
  if (is.infinite(log_med_ratio)) {
    log_med_ratio <- 0
  }
  log_mean_ratio <- log2(mean(abund_1) / mean(abund_2))
  if (is.nan(log_mean_ratio)) {
    log_mean_ratio <- 0
  }
  if (is.infinite(log_mean_ratio)) {
    log_mean_ratio <- 0
  }
  list(
    log2_median_ratio = log_med_ratio,
    log2_mean_ratio = log_mean_ratio,
    median_diff = median(abund_1) - median(abund_2),
    mean_diff = mean(abund_1) - mean(abund_2),
    mean_treat1 = mean(abund_1),
    mean_treat2 = mean(abund_2),
    wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value,
    hartigan_dip_treat1 = diptest::dip.test(abund_1)$p.value,
    hartigan_dip_treat2 = diptest::dip.test(abund_2)$p.value,
    bimodality_coeff_treat1 = modes::bimodality_coefficient(abund_1),
    bimodality_coeff_treat2 = modes::bimodality_coefficient(abund_2)
  )
}

initial_filter <- function(obj, ambiguous_filter, kingdom_filter) {
  # Check to see if filtering is necessary and then return the metacoder/taxmap object
  # TODO: Make this function more mature.
  # Filter ambiguous taxon names, and filter kingdoms and return unfiltered taxmap
  if (!is.null(ambiguous_filter) & !is.null(kingdom_filter)) {
    filtered_meta <- obj
    for (term in ambiguous_filter) {
      filtered_meta <- filtered_meta %>%
        taxa::filter_taxa(!stringr::str_detect(taxon_names, term))
    }
    for (term in kingdom_filter) {
      filtered_meta <- filtered_meta %>%
        taxa::filter_taxa(taxon_names == term, subtaxa = TRUE)
    }
    return(filtered_meta)

    # Filter ambiguous taxon names and return the filtered taxmap
  } else if (!is.null(ambiguous_filter) & is.null(kingdom_filter)) {
    filtered_meta <- obj
    for (term in ambiguous_filter) {
      filtered_meta <- filtered_meta %>%
        taxa::filter_taxa(!stringr::str_detect(taxon_names, term))
    }
    return(filtered_meta)

    # Filter kingdoms and return the filtered taxmap
  } else if (!is.null(kingdom_filter) & is.null(ambiguous_filter)) {
    filtered_meta <- obj
    for (term in kingdom_filter) {
      filtered_meta <- filtered_meta %>%
        taxa::filter_taxa(taxon_names == term, subtaxa = TRUE)
    }
    return(filtered_meta)

    # If filtering is unceccessary return unfiltered taxmap
  } else {
    return(obj)
  }
}
