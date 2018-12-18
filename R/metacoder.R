#' @title Filter Samples from Metacoder Objects
#' @description This function provides a flexible way to filter unwanted samples from the
#' \emph{"otu_abundance"} and \emph{"sample_data"} observations of a MicrobiomeR formatted object.
#' @param obj A Taxmap/metacoder object.
#' @param .f_transform A function used for transforming the data.  Default: NULL
#' @param .f_filter A function used for summarising the data like 'sum' or 'mean'.  Default: NULL
#' @param .f_condition A function that takes the summarised data and applied a condition like x > 10000.  Default: NULL
#' @param validated This parameter provides a way to override validation steps.  Use carefully.  Default: FALSE
#' @param ... An optional list of parameters to use in the .f_filter function specified
#' @return Returns a metacoder object with samples that pass the filters.
#' @details Get the samples to keep by using purr and the user supplied transform and filter + condition formulas.
#' The purr packge allows the use of anonymous functions as described in the link below:
#'
#' \url{https://jennybc.github.io/purrr-tutorial/ls03_map-function-syntax.html#anonymous_function,_formula}
#' @examples
#' \dontrun{
#' if(interactive()){
#' # Use the sample_filter early on in your analysis
#' library(MicrobiomeR)
#' library(metacoder)
#' library(taxa)
#' # Convert Phyloseq object to metacoder object
#' metacoder_obj <- as_MicrobiomeR_format(obj = phyloseq_obj, format = "raw_format")
#'
#' # Remove Archaea from the metacoder object
#' metacoder_obj <- filter_taxa(
#'   obj = metacoder_obj,
#'   taxon_names == "Archaea",
#'   subtaxa = TRUE,
#'   invert = TRUE)
#'
#' # Ambiguous Annotation Filter - Remove taxonomies with ambiguous names
#' metacoder_obj <- filter_ambiguous_taxa(metacoder_obj, subtaxa = TRUE)
#'
#' # Low Sample Filter - Remove the low samples
#' metacoder_obj <- sample_filter(obj          = metacoder_obj,
#'                                .f_filter    = ~sum(.),
#'                                .f_condition = ~.>= 20, validated = TRUE)
#'  }
#' }
#' @export
#' @family Metacoder Filters
#' @rdname sample_filter
#' @seealso
#'  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}},  \code{\link[MicrobiomeR]{transformer}}
#'
#'  \code{\link[taxa]{get_dataset}}
#'
#'  \code{\link[dplyr]{select_all}},  \code{\link[dplyr]{select}},  \code{\link[dplyr]{filter}}
#'
#'  \code{\link[purrr]{map}},  \code{\link[purrr]{keep}}
#' @importFrom taxa get_dataset
#' @importFrom dplyr select_if select filter
#' @importFrom purrr map keep
sample_filter <- function(obj, .f_transform = NULL, .f_filter = NULL, .f_condition = NULL, validated = FALSE, ...) {
  mo_clone <- obj$clone()
  # Make sure the required functions are provided.
  if (!is.null(.f_filter) && !is.null(.f_condition)) {

    # Validate and get otu_abundance data
    mo_clone <- validate_MicrobiomeR_format(obj = mo_clone, valid_formats = c("raw_format", "basic_format"),
                                            force_format = TRUE, validated = validated, min_or_max = min)
    abund_data <- taxa::get_dataset(obj = mo_clone, data = "otu_abundance")

    #  Get the sample data columns to work with
    if (is.null(.f_transform)) { # Get raw sample data
      sample_cols <- abund_data %>% dplyr::select_if(is.numeric)
    } else { # Get transfmormed sample data
      sample_cols <- transformer(.data = abund_data, func = .f_transform) %>% dplyr::select_if(is.numeric)
    }
    # Get the samples to keep by using purr and the user supplied filter and condition formulas
    samples_to_keep <- purrr::map(sample_cols, .f_filter, list(...)) %>% # Apply a summary function like 'sum' or 'mean'
      purrr::map(.f_condition) %>% # Apply a function to the sample data that does some comparison (f(x)x>10000)
      purrr::keep(~ . == TRUE) %>%
      names() # Determine which samples to keep

    # Update the otu_abundance and sample_data in the metacoder object by removing samples
    other_vars <- abund_data %>% dplyr::select_if(function(x) is.numeric(x) == FALSE) %>% colnames() #%>% purrr::discard(~.=="taxon_id")
    mo_clone$data$otu_abundance  <- dplyr::select(abund_data, c(other_vars, samples_to_keep))
    mo_clone$data$sample_data <- mo_clone$data$sample_data %>% dplyr::filter(sample_id %in% samples_to_keep)
    return(mo_clone)
  } else {
    stop("You have to supply a filter formula AND a condition formula.")
  }
}

#' @title Filter Taxon Ids from Metacoder Objects
#' @description This function provides a flexible way to filter unwanted taxon_ids from the taxmap object and from the
#' observations of a MicrobiomeR formatted object.
#' @inheritParams sample_filter
#' @return Returns a metacoder object with taxon_ids that pass the filters.
#' @pretty_print TRUE
#' @details Get the taxon_ids to keep by using purr and the user supplied transform and filter + condition formulas.
#' The purr packge allows the use of anonymous functions as described in the link below:
#'
#' \url{https://jennybc.github.io/purrr-tutorial/ls03_map-function-syntax.html#anonymous_function,_formula}
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(MicrobiomeR)
#' library(metacoder)
#' library(taxa)
#' # Convert Phyloseq object to metacoder object
#' metacoder_obj <- as_MicrobiomeR_format(obj = phyloseq_obj, format = "raw_format")
#'
#' # Remove Archaea from the metacoder object
#' metacoder_obj <- filter_taxa(
#'   obj = metacoder_obj,
#'   taxon_names == "Archaea",
#'   subtaxa = TRUE,
#'   invert = TRUE)
#'
#' # Ambiguous Annotation Filter - Remove taxonomies with ambiguous names
#' metacoder_obj <- filter_ambiguous_taxa(metacoder_obj, subtaxa = TRUE)
#'
#' # Low taxa filter - Remove the low taxon_ids
#' metacoder_obj <- taxon_id_filter(obj        = metacoder_obj,
#'                                .f_filter    = ~sum(.),
#'                                .f_condition = ~.>= 2000, validated = TRUE)
#'  }
#' }
#' @export
#' @family Metacoder Filters
#' @rdname taxon_id_filter
#' @seealso
#'  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}},  \code{\link[MicrobiomeR]{transposer}},  \code{\link[MicrobiomeR]{transformer}}
#'
#'  \code{\link[taxa]{get_dataset}},  \code{\link[taxa]{filter_taxa}}
#'
#'  \code{\link[dplyr]{select_all}}
#'
#'  \code{\link[purrr]{map}},  \code{\link[purrr]{keep}}
#' @importFrom taxa get_dataset filter_taxa
#' @importFrom dplyr select_if
#' @importFrom purrr map keep
taxon_id_filter <- function(obj, .f_transform = NULL, .f_filter = NULL, .f_condition = NULL, validated = FALSE, ...) {
  mo_clone <- obj$clone()
  # Make sure the required functions are provided.
  if (!is.null(.f_filter) && !is.null(.f_condition)) {

    # Validate and get taxa_abundance data
    mo_clone <- validate_MicrobiomeR_format(obj = mo_clone, valid_formats = c("raw_format", "basic_format"),
                                            force_format = TRUE, validated = validated, min_or_max = min)
    abund_data <- taxa::get_dataset(obj = mo_clone, data = "taxa_abundance")
    #  Get the taxon_id data columns to work with
    if (is.null(.f_transform)) { # Get raw sample data
      taxon_id_cols <- transposer(abund_data, ids = "taxon_id", header_name = "samples", preserved_categories = FALSE) %>%
        dplyr::select_if(is.numeric)
    } else { # Get transfmormed sample data
      taxon_id_cols <- transformer(.data = abund_data, func = .f_transform) %>%
        transposer(ids = "taxon_id", header_name = "samples", preserved_categories = FALSE) %>%
        dplyr::select_if(is.numeric)
    }
    # Get the samples to keep by using purr and the user supplied filter and condition formulas
    # The purr packge allows the use of anonymous functions as described in the link below:
    # https://jennybc.github.io/purrr-tutorial/ls03_map-function-syntax.html#anonymous_function,_formula
    taxon_ids_to_keep <- purrr::map(taxon_id_cols, .f_filter, ...) %>% # Apply a summary function like 'sum' or 'mean'
      purrr::map(.f_condition) %>% # Apply a function to the sample data that does some comparison (f(x)x>10000)
      purrr::keep(~ . == TRUE) %>%
      names() # Determine which samples to keep
    # Update the observation data tables and the taxmap object
    mo_clone <- taxa::filter_taxa(mo_clone, taxon_ids %in% taxon_ids_to_keep, reassign_obs = FALSE)
    return(mo_clone)
  } else {
    stop("You have to supply a filter formula AND a condition formula.")
  }
}

#' @title Filter OTU Ids from Metacoder Objects
#' @description This function provides a flexible way to filter unwanted otu_ids from the taxmap object and from the
#' observations of a MicrobiomeR formatted object.
#' @inheritParams sample_filter
#' @return Returns a metacoder object with otu_ids that pass the filters.
#' @pretty_print TRUE
#' @details Get the otu_ids to keep by using purr and the user supplied transform and filter + condition formulas.
#' The purr packge allows the use of anonymous functions as described in the link below:
#'
#' \url{https://jennybc.github.io/purrr-tutorial/ls03_map-function-syntax.html#anonymous_function,_formula}
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Metacoder Filters
#' @rdname otu_id_filter
#' @seealso
#'  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}},  \code{\link[MicrobiomeR]{transposer}},  \code{\link[MicrobiomeR]{transformer}}
#'
#'  \code{\link[taxa]{get_dataset}},  \code{\link[taxa]{filter_taxa}}
#'
#'  \code{\link[dplyr]{select_all}}
#'
#'  \code{\link[purrr]{map}},  \code{\link[purrr]{keep}}
#' @importFrom taxa get_dataset filter_taxa
#' @importFrom dplyr select_if
#' @importFrom purrr map keep
otu_id_filter <- function(obj, .f_transform = NULL, .f_filter = NULL, .f_condition = NULL, validated = FALSE, ...) {
  mo_clone <- obj$clone()
  # Make sure the required functions are provided.
  if (!is.null(.f_filter) && !is.null(.f_condition)) {
    # Validate and get taxa_abundance data
    mo_clone <- validate_MicrobiomeR_format(obj = mo_clone, valid_formats = c("raw_format", "basic_format"),
                                                         force_format = TRUE, validated = validated, min_or_max = min)
    fmt <- which_format(mo_clone)
    abund_data <- taxa::get_dataset(obj = mo_clone, data = "otu_abundance")
    #  Get the taxon_id data columns to work with
    if (is.null(.f_transform)) { # Get raw sample data
      otu_id_cols <- transposer(abund_data, ids = "otu_id", header_name = "samples", preserved_categories = FALSE) %>%
        dplyr::select_if(is.numeric)
    } else { # Get transfmormed sample data
      otu_id_cols <- transformer(.data = abund_data, func = .f_transform) %>%
        transposer(ids = "otu_id", header_name = "samples", preserved_categories = FALSE) %>%
        dplyr::select_if(is.numeric)
    }
    # Get the samples to keep by using purr and the user supplied filter and condition formulas
    otu_ids_to_keep <- purrr::map(otu_id_cols, .f_filter, ...) %>% # Apply a summary function like 'sum' or 'mean'
      purrr::map(.f_condition) %>% # Apply a function to the sample data that does some comparison (f(x)x>10000)
      purrr::keep(~ . == TRUE) %>%
      names() # Determine which samples to keep
    # Get an updated vector of the observation tables to filter
    otu_table_list <- pkg.private$format_table_list$otu_tables
    otu_table_list <- otu_table_list[otu_table_list %in% names(mo_clone$data)]
    # Update the observation data tables and the taxmap object
    mo_clone <- taxa::filter_obs(mo_clone, otu_table_list, otu_id %in% otu_ids_to_keep, drop_taxa = TRUE)
    return(mo_clone)
  } else {
    stop("You have to supply a filter formula AND a condition formula.")
  }
}

#' @title Agglomerate Metacoder Objects
#' @description A function similar to the \code{\link[phyloseq:tax_glom]{phyloseq::tax_glom}} function,
#' that assembles abundance data at a specified rank.  This removes subtaxa and reassigns the
#' values at the specified rank.
#' @param obj A Taxmap/metacoder object.
#' @param rank The rank that will be agllomerated to.
#' @param validated This parameter provides a way to override validation steps.  Use carefully.  Default: FALSE
#' @return A taxmap object that has been agglomerated at the specified rank.
#' @pretty_print TRUE
#' @details This function helps analyzing taxonomic data at intermediate ranks by agglomerating the observation data
#' and the taxmap object.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Metacoder Filters
#' @rdname agglomerate_metacoder
#' @seealso
#'  \code{\link[taxa]{filter_taxa}}
#'
#'  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}}
#' @importFrom taxa filter_taxa taxon_ranks
agglomerate_metacoder <- function(obj, rank, validated = FALSE) {
  mo_clone <- obj$clone()
  mo_clone <- validate_MicrobiomeR_format(obj = mo_clone, valid_formats = c("raw_format", "basic_format"),
                                          force_format = TRUE, validated = validated, min_or_max = min)
  # Agglomerate
  mo_clone <- taxa::filter_taxa(mo_clone, taxon_ranks == rank,
                                supertaxa = TRUE, reassign_obs = FALSE
  )
  return(mo_clone)
}

#' @title OTU Proportion FIlter
#' @description This function filters OTU values from the observation data and the taxmap object
#' based on a minimum proportional mean across samples per OTU.
#' @param obj A Taxmap/metacoder object.
#' @param otu_percentage The minimum percentage used to compare against the proportional OTU mean.  Default: 5e-05
#' @param validated This parameter provides a way to override validation steps.  Use carefully.  Default: FALSE
#' @return Returns a taxmap object that contains otu_ids that have passed the above filter.
#' @details This type of filtering is used to remove OTUs that do not have a specified mean proportion.
#' This function must be used conservatively, hence why the default otu_percentage is so low.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Metacoder Filters
#' @rdname otu_proportion_filter
#' @seealso
#'  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}},\code{\link[MicrobiomeR]{otu_id_filter}}
otu_proportion_filter <- function(obj, otu_percentage = 0.00005, validated = FALSE) {
  mo_clone <- obj$clone()
  mo_clone <- validate_MicrobiomeR_format(obj = mo_clone, valid_formats = c("raw_format", "basic_format"),
                                                       force_format = TRUE, validated = validated, min_or_max = min)
  mo_clone <- otu_id_filter(obj = mo_clone, .f_transform = ~./sum(.), .f_filter = ~mean(.), .f_condition = ~.> otu_percentage)
  return(mo_clone)
}

#' @title OTU Prevalence Filter (Metacoder)
#' @description This function filters observations by thier prevelance across samples.
#' @param obj A Taxmap/metacoder object.
#' @param minimum_abundance The minimum abundance needed per observation per sample.  Default: 5
#' @param rel_sample_percentage The percentage of samples per observation that meet the minimum abundance.  Default: 0.5
#' @param validated This parameter provides a way to override validation steps.  Use carefully.  Default: FALSE
#' @return Returns a taxmap object that contains taxon_ids that have passed the above filter.
#' @pretty_print TRUE
#' @details The otu_prevalence_filter filters taxon_ids that do not appear more than a certian amount of times (minimum abundance) in a certain percentage of
#' samples (rel_sample_percentage). The \href{http://web.stanford.edu/class/bios221/MicrobiomeWorkflowII.html#filtering}{phyloseq workflow} calls for a minimum abundance of 5 across %50 of the samples.
#' This filtering method is considered unsupervised, because it soley relies on the data in this experiment (OTU ids).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Metacoder Filters
#' @rdname otu_prevalence_filter
#' @seealso
#'  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}}
#'  \code{\link[metacoder]{calc_prop_samples}}
#'  \code{\link[dplyr]{filter}}
#'  \code{\link[taxa]{filter_taxa}}
#' @importFrom metacoder calc_prop_samples
#' @importFrom dplyr filter
#' @importFrom taxa filter_taxa
otu_prevalence_filter <- function(obj, minimum_abundance = 5, rel_sample_percentage = 0.5,
                                  validated = FALSE) {
  mo_clone <- obj$clone()
  mo_clone <- validate_MicrobiomeR_format(obj = mo_clone, valid_formats = c("basic_format"),
                                          force_format = TRUE, validated = validated, min_or_max = min)
  # Calculate the ids that need to be removed
  ids_to_remove <- metacoder::calc_prop_samples(mo_clone, "taxa_abundance", more_than = minimum_abundance) %>% # Calculate sample proportions per taxa with min abundance
    dplyr::filter(n_samples < rel_sample_percentage) # Filter samples with less than the relative sample percentage
  # Prevalence Filtering
  mo_clone <- taxa::filter_taxa(mo_clone, !taxon_ids %in% ids_to_remove$taxon_id, reassign_obs = FALSE)
  return(mo_clone)
}

#' @title Taxonomic Prevalence Filter (Metacoder)
#' @description This function filters observations at a specific rank by thier prevelance across samples.
#' @param obj A Taxmap/metacoder object.
#' @param rank The rank being analyzed for prevalence across samples.
#' @param minimum_abundance The minimum abundance needed per observation per sample.  Default: 5
#' @param rel_sample_percentage The percentage of samples per observation that meet the minimum abundance.  Default: 0.5
#' @param validated This parameter provides a way to override validation steps.  Use carefully.  Default: FALSE
#' @return Returns a taxmap object that contains taxon_ids that have passed the above fiter.
#' @pretty_print TRUE
#' @details The taxa_prevalence_filter filters taxon_ids that do not appear more than a certian amount of times (minimum abundance) in a certain percentage of
#' samples (rel_sample_percentage) at the specified agglomerated rank (rank). The \href{http://web.stanford.edu/class/bios221/MicrobiomeWorkflowII.html#filtering}{phyloseq workflow} calls for a minimum abundance of 5 across %50 of the samples.
#' This method is considered supervised, because the filtering is done based on taxonomic annotation (taxon_ids), which is assigned based on a reference database (SILVA or GreenGenes).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Metacoder Filters
#' @rdname taxa_prevalence_filter
#' @seealso
#'  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}}
#'  \code{\link[metacoder]{calc_prop_samples}}
#'  \code{\link[dplyr]{filter}}
#'  \code{\link[taxa]{filter_taxa}}
#' @importFrom metacoder calc_prop_samples
#' @importFrom dplyr filter
#' @importFrom taxa filter_taxa
taxa_prevalence_filter <- function(obj, rank, minimum_abundance = 5, rel_sample_percentage = 0.5,
                                   validated = FALSE) {
  mo_clone <- obj$clone()
  mo_clone <- validate_MicrobiomeR_format(obj = mo_clone, force_format = TRUE, validated = validated,
                                          min_or_max = min, valid_formats = c("basic_format"))
  # Calculate the ids that need to be removed based on taxonomic rank
  ids_to_remove <- agglomerate_metacoder(obj = mo_clone, rank = rank, validated = TRUE) %>% # Agglomeration
    metacoder::calc_prop_samples("taxa_abundance", more_than = minimum_abundance) %>% # Calculate sample proportions per taxa with min abundance
    dplyr::filter(n_samples < rel_sample_percentage) # Filter samples with less than the relative sample percentage
  # Taxonomic Prevalence Filtering
  mo_clone <- taxa::filter_taxa(mo_clone, !taxon_ids %in% ids_to_remove$taxon_id, reassign_obs = FALSE)
  return(mo_clone)
}

#' @title Coefficient of Variation Filter
#' @description This function filters OTUs that have a variance higher than the
#' specified CoV.
#' @param obj A Taxmap/metacoder object.
#' @param coefficient_of_variation The maximum CoV that an OTU can have.
#' @param validated This parameter provides a way to override validation steps.  Use carefully.  Default: FALSE
#' @return Returns a taxmap object that contains otu_ids that have passed the above filter.
#' @pretty_print TRUE
#' @details This function helps remove OTUs that have an unusually high variance using the coefficient
#' of variation.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Metacoder Filters
#' @rdname cov_filter
#' @seealso
#'  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}},  \code{\link[MicrobiomeR]{otu_id_filter}}
#'  \code{\link[dplyr:summarise_all]{summarise_if}}
#' @importFrom dplyr summarise_if
#' @importFrom stats sd
cov_filter <- function(obj, coefficient_of_variation, validated = FALSE) {
  mo_clone <- obj$clone()
  mo_clone <- validate_MicrobiomeR_format(obj = mo_clone, valid_formats = c("raw_format", "basic_format"),
                                          force_format = TRUE, validated = validated, min_or_max = min)
  # Standardize abundances to the median sequencing depth
  total <- mo_clone$data$otu_abundance %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    as.numeric() %>%
    median()
  standf <- function(x, t = total) round(t * (x / sum(x)))
  # Filter OTUs that don't pass the maximum coefficient of variation.
  mo_clone <- otu_id_filter(obj = mo_clone, .f_transform = standf, .f_filter = ~sd(.)/mean(.), .f_condition = ~.<coefficient_of_variation)
  return(mo_clone)
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
#'
#'  \code{\link[modes]{bimodality_coefficient}}
#' @importFrom diptest dip.test
#' @importFrom modes bimodality_coefficient
#' @importFrom stats wilcox.test median
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
