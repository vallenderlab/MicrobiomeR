


#' @title Get Metacoder Object
#' @description Create a metacoder object using a phyloseq object or rdata file.
#' @param phyloseq_object A phyloseq object.  Default: NULL
#' @param rdata_file A .Rdata file.  Default: NULL
#' @param filter_flag A flag that will determine if the metacoder object gets filtered.  Default: FALSE
#' @param ambiguous_filter A vector denoting ambiguous taxonomies to filter.  Default: c("uncultured", "Unassigned", "Ambiguous")
#' @param kingdom_filter A character vector denoting which Kingodm's to keep. Default: c("Bacteria")
#' @param workflow_func A function that takes a metacoder object and a comparison function (See \code{\link{metacoder_workflow_1}})
#' @param comp_func A function used for metacoders compare_groups function.  (See \code{\link{metacoder_comp_func_1}})
#' @return Returns a metacoder object to be used for generating a heat tree.
#' @details Creates a metacoder object that by default compares treatment groups for downstream analysis.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Data Importer
#' @rdname get_metacoder_obj
#' @seealso
#'  \code{\link[metacoder]{parse_phyloseq}},\code{\link[metacoder]{calc_obs_props}},\code{\link[metacoder]{calc_taxon_abund}},\code{\link[metacoder]{compare_groups}}
#' @importFrom metacoder parse_phyloseq calc_obs_props calc_taxon_abund compare_groups
get_metacoder_obj <- function(phyloseq_object=NULL, rdata_file=NULL, filter_flag = FALSE,
                              ambiguous_filter=c("uncultured", "Unassigned", "Ambiguous"),
                              kingdom_filter=c("Bacteria"), workflow_func = metacoder_workflow_1,
                              comp_func = metacoder_comp_func_1) {
  if (!is.null(rdata_file)){
    # Load data
    load(rdata_file)
    # Reassign loaded metacoder object to "master_metacoder"
    metacoder_obj <- get("metacoder_obj")
  } else {
    # Convert Phyloseq objects to Metacoder objects using the workflow/comparison funcctions
    metacoder_obj <- metacoder::parse_phyloseq(phyloseq_object)
    if (!is.null(workflow_func)){
      metacoder_obj <- workflow_func(metacoder_obj = metacoder_obj, func = comp_func)
    }
  }

  if (filter_flag == TRUE) {
    # Put data through initial filtering of ambiguous and kingdom data
    # Note always do this after calculating relative/proportionate abundance
    processed_metacoder <- metacoder_initial_filter(
      metacoder_obj = metacoder_obj,
      ambiguous_filter = ambiguous_filter,
      kingdom_filter = kingdom_filter
    )
    return(processed_metacoder)
  } else {
    return(metacoder_obj)
  }
}

#' @title Metacoder Workflow Function #1
#' @description A function for analyzing or manipulating metacoder taxmap objects.
#' @param metacoder_obj A taxmap generated by metacoder.
#' @param func A comparison function for use with \code{\link[metacoder]{compare_groups}} Default: metacoder_comp_func_1
#' @return Returns a altered metacoder object.
#' @pretty_print TRUE
#' @details Takes a metacoder object and a function for manipulating metacoder objects.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Data Manipulators
#' @rdname metacoder_workflow_1
#' @seealso
#'  \code{\link[metacoder]{calc_obs_props}},\code{\link[metacoder]{calc_taxon_abund}},\code{\link[metacoder]{compare_groups}}
#' @importFrom metacoder calc_obs_props calc_taxon_abund compare_groups
metacoder_workflow_1 <- function(metacoder_obj, func=metacoder_comp_func_1) {
  # Calculate a proportionate abundance value for each sample
  metacoder_obj$data$otu_table <- metacoder::calc_obs_props(metacoder_obj, data = "otu_table", cols = metacoder_obj$data$sample_data$sample_id, other_cols = TRUE)
  # Calculate the proportionate abundace per-taxon value for each sample
  metacoder_obj$data$tax_table <- metacoder::calc_taxon_abund(metacoder_obj, data = "otu_table", cols = metacoder_obj$data$sample_data$sample_id)
  # Calculate the differences between Stress vs Control
  metacoder_obj$data$diff_table <- metacoder::compare_groups(metacoder_obj, data = "tax_table", cols = metacoder_obj$data$sample_data$sample_id, groups = metacoder_obj$data$sample_data[["X.TreatmentGroup"]], func = func, other_cols = TRUE)
  # Rename the tax_table columns for future filtering or manipulations
  # TODO: Find out if the following code is still relevant or can be updated.
  colnames(metacoder_obj$data$tax_table) <- sapply(colnames(metacoder_obj$data$tax_table), function(x) if(x!="taxon_id"){sprintf("T%s", x)}else{x})
}

#' @title Metacoder Comparison Function #1
#' @description A comparison function for metacoder::compare_groups "func" parameter.
#' @param abund_1 A character vector of abundances.
#' @param abund_2 A character vector of abundances.
#' @return A list of statistical results used to compare groups.
#' @pretty_print TRUE
#' @details This function is used by metacoder::compare_groups in order to compare
#' every combination of treatment groups.
#' @export
#' @family Data Manipulators
#' @rdname metacoder_comp_func_1
#' @seealso
#'  \code{\link[diptest]{dip.test}}
#'  \code{\link[modes]{bimodality_coefficient}}
#' @importFrom diptest dip.test
#' @importFrom modes bimodality_coefficient
metacoder_comp_func_1 <- function(abund_1, abund_2) {
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

#' @title Metacoder Inital Filter
#' @description A function for filtering ambiguous taxonomy annotations and unwanted Kindgom level taxonomies.
#' @param metacoder_obj A metacoder object.
#' @param ambiguous_filter A vector of ambiguous taxonomy annotations to filter out.
#' @param kingdom_filter A vector of Kingodm level taxonomies to keep.
#' @return Returns a filtered metacoder object.
#' @pretty_print TRUE
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Filters
#' @rdname metacoder_initial_filter
#' @seealso
#'  \code{\link[taxa]{filter_taxa}}
#'  \code{\link[stringr]{str_detect}}
#' @importFrom taxa filter_taxa
#' @importFrom stringr str_detect
metacoder_initial_filter <- function(metacoder_obj, ambiguous_filter, kingdom_filter) {
  # Check to see if filtering is necessary and then return the metacoder/taxmap object
  # TODO: Make this function more mature.
  # Filter ambiguous taxon names, and filter kingdoms and return unfiltered taxmap
  if (!is.null(ambiguous_filter) & !is.null(kingdom_filter)) {
    filtered_meta <- metacoder_obj
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
    filtered_meta <- metacoder_obj
    for (term in ambiguous_filter) {
      filtered_meta <- filtered_meta %>%
        taxa::filter_taxa(!stringr::str_detect(taxon_names, term))
    }
    return(filtered_meta)

    # Filter kingdoms and return the filtered taxmap
  } else if (!is.null(kingdom_filter) & is.null(ambiguous_filter)) {
    filtered_meta <- metacoder_obj
    for (term in kingdom_filter) {
      filtered_meta <- filtered_meta %>%
        taxa::filter_taxa(taxon_names == term, subtaxa = TRUE)
    }
    return(filtered_meta)

    # If filtering is unceccessary return unfiltered taxmap
  } else {
    return(metacoder_obj)
  }
}