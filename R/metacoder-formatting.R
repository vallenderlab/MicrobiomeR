#' @title Which MicrobiomeR Format
#' @description A function for looking at a metacoder object and returning the identified MicrobiomeR format.
#' @param obj A Taxmap/metacoder object.
#' @return If the format is verified it returns a character string denoting the identified format.
#' @pretty_print TRUE
#' @details This function is used to get basic information about the format of the taxmap object
#' that is supplied by the user.
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(MicrobiomeR)
#' fmt <- which fmt(MicrobiomeR::raw_silva)
#' print(fmt)
#'  }
#' }
#' @export
#' @family Validation
#' @rdname which_format
#' @importFrom crayon yellow red
which_format <- function(obj) {
  mo_clone <- obj$clone()
  # Table names
  table_names <- names(mo_clone$data)
  # Format names
  raw_names <- pkg.private$format_table_list$raw_format
  basic_names <- pkg.private$format_table_list$basic_format
  analyzed_names <- pkg.private$format_table_list$analyzed_format
  phyloseq_names <- pkg.private$format_table_list$phyloseq_format
  # Format Flags
  raw_flag <- all(raw_names %in% table_names)
  basic_flag <- all(basic_names %in% table_names)
  analyzed_flag <- all(analyzed_names %in% table_names)
  phyloseq_flag <- all(phyloseq_names %in% table_names)
  if (phyloseq_flag) {
    other_flag <- any(c(raw_flag, basic_flag, analyzed_flag))
    if (!other_flag) {
      message(crayon::yellow("Your object is in the phyloseq format!"))
      message(crayon::yellow("Please format your metacoder object to continue analysis."))
      return("phyloseq_format")
    } else {
      message(crayon::yellow(sprintf("The table names in the metacoder object are: %s", paste(table_names, collapse = ", "))))
      message(crayon::yellow("You have a mix between phyloseq format and other format."))
      return("mixed_format")
    }
  } else
    # Logic for returning the format type.
    if (analyzed_flag){
      return("analyzed_format")
    } else if (basic_flag) {
      return("basic_format")
    } else if (raw_flag) {
      return("raw_format")
    } else  {
      warning(crayon::red(sprintf("The table names in the metacoder object are: %s", paste(table_names, collapse = ", "))))
      warning(crayon::red("The object is not in a recognized format."))
      return("unknown_format")
    }
}

#' @title Is Raw MicrobiomeR Format
#' @description This function returns a logical based on weather or not the object is in the raw_format.
#' @param obj A Taxmap/metacoder object.
#' @return A logical (TRUE/FALSE).
#' @pretty_print TRUE
#' @details The "raw_format" is Level 1. in the [MicrobiomeR_Formats] hierarchy.
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(MicrobiomeR)
#' if(is_raw_format(MicrobiomeR::raw_silva)){
#'     print("It's raw!")
#' } else {
#'     fmt <- which_format(MicrobiomeR::raw_silva)
#'     print("It's not analyzed!")
#'     print(sprintf("It's %s!", fmt))
#' }
#'  }
#' }
#' @export
#' @family Validation
#' @rdname is_raw_format
is_raw_format <- function(obj) {
  fmt <- which_format(obj)
  if (fmt == "raw_format"){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @title Is Basic MicrobiomeR Format
#' @description This function returns a logical based on weather or not the object is in the basic_format.
#' @param obj A Taxmap/metacoder object.
#' @return A logical (TRUE/FALSE).
#' @pretty_print TRUE
#' @details The "basic_format" is Level 2. in the [MicrobiomeR_Formats] hierarchy.
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(MicrobiomeR)
#' if(is_basic_format(MicrobiomeR::raw_silva)){
#'     print("It's basic!")
#' } else {
#'     fmt <- which_format(MicrobiomeR::raw_silva)
#'     print("It's not analyzed!")
#'     print(sprintf("It's %s!", fmt))
#' }
#'  }
#'}
#' @export
#' @family Validation
#' @rdname is_basic_format
is_basic_format <- function(obj) {
  fmt <- which_format(obj)
  if (fmt == "basic_format"){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @title Is Analyzed MicrobiomeR Format
#' @description This function returns a logical based on weather or not the object is in the analyzed_format.
#' @param obj A Taxmap/metacoder object.
#' @return A logical (TRUE/FALSE).
#' @pretty_print TRUE
#' @details The "analyzed_format" is Level 3. in the [MicrobiomeR_Formats] hierarchy.
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(MicrobiomeR)
#' if(is_analyzed_format(MicrobiomeR::raw_silva)){
#'     print("It's analyzed!")
#' } else {
#'     fmt <- which_format(MicrobiomeR::raw_silva)
#'     print("It's not analyzed!")
#'     print(sprintf("It's %s!", fmt))
#' }
#'  }
#' }
#' @export
#' @family Validation
#' @rdname is_analyzed_format
is_analyzed_format <- function(obj) {
  fmt <- which_format(obj)
  if (fmt == "analyzed_format"){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' @title Is Phyloseq MicrobiomeR Format
#' @description This function returns a logical based on weather or not the object is in the phyloseq_format.
#' @param obj A Taxmap/metacoder object.
#' @return A logical (TRUE/FALSE).
#' @pretty_print TRUE
#' @details The "phyloseq_format" is Level 0. in the [MicrobiomeR_Formats] hierarchy.
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(MicrobiomeR)
#' if(is_analyzed_format(MicrobiomeR::raw_silva)){
#'     print("It's phyloseq format!")
#' } else {
#'     fmt <- which_format(MicrobiomeR::raw_silva)
#'     print("It's not phyloseq fomrat!")
#'     print(sprintf("It's %s!", fmt))
#' }
#'  }
#' }
#' @export
#' @family Validation
#' @rdname is_phyloseq_format
is_phyloseq_format <- function(obj) {
  fmt <- which_format(obj)
  if (fmt == "phyloseq_format"){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @title Order Metacoder Observation Data
#' @description A function for changing the order of the observation data in a metacoder object
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @return A taxmap object with observation data in the proper order for downstream analysis.
#' @pretty_print TRUE
#' @details Changes the order of the observation tables in \strong{metacoder_object$data} to
#' otu_abundance, otu_annotations, otu_proportions, sample_data, phy_tree, taxa_abundance,
#' taxa_proportions, statistical_data, and stats_tax_data respectively.
#' @export
#' @family Validation
#' @rdname order_metacoder_data
#' @importFrom crayon silver
order_metacoder_data <- function(obj) {
  obj <- object_handler(obj = obj)
  mo_clone <- obj$clone()
  # Create list of expected table names
  expected_names <- pkg.private$format_table_list$expected_table_order

  # Create vector that contains expected tables names already in the metacoder object
  table_order <- names(expected_names[names(expected_names) %in% names(mo_clone$data)])
  # Get the table names that aren't expected but in the metacoder object
  other_names <- names(mo_clone$data)[!names(mo_clone$data) %in% c(expected_names)]
  table_order <- c(table_order, other_names)
  mo_clone$data <- mo_clone$data[table_order]
  message(crayon::silver("Ordered Data."))
  return(mo_clone)
}

#' @title Validate MicrobiomeR Format
#' @description This funciton validates that taxmap/metacoder objects are in a valid format MicrobiomeR format.
#' @param obj A Taxmap/metacoder object.
#' @param validated This parameter provides a way to override validation steps.  Use carefully.  Default: FALSE
#' @param valid_formats A vector of formats that are used for validation.
#' @param force_format A logical denoting if the selected format is to be forced.  Default: FALSE
#' @param min_or_max A base function (min or max) that determines which format is to be forced.
#' This is particularly useful if you provide multiple \emph{valid_formats}.  Min will choose the lowest
#' level format, while max will choose the highest level format.  Default: base::max
#' @param ... An optional list of parameters to use in \code{\link[MicrobiomeR]{as_MicrobiomeR_format}}.
#' @return If the object is validated, a Taxmap/metacoder object.
#' @details This function can provide a way to check if a taxmap object has undergone a
#' \code{\link[MicrobiomeR:MicrobiomeR_Workflow]{MicrobiomeR Style Workflow}}.
#' @export
#' @family Validation
#' @rdname validate_MicrobiomeR_format
#' @seealso
#'  \code{\link[MicrobiomeR]{which_format}},  \code{\link[MicrobiomeR]{as_MicrobiomeR_format}}
#' @importFrom glue glue
#' @importFrom crayon bgWhite green yellow red
validate_MicrobiomeR_format <- function(obj, validated = FALSE, valid_formats, force_format = FALSE, min_or_max = base::max, ...) {
  mo_clone <- obj$clone()
  format_list <- pkg.private$format_level_list
  if (validated == TRUE) {
    return(mo_clone)
  }
  fmt <- which_format(obj = mo_clone)
  rank_list <- c(format_list[[fmt]])
  high_rank <- fmt
  if (fmt %in% valid_formats) {
    return(mo_clone)
  } else if (force_format == TRUE) {
    for (v_fmt in valid_formats) {
      rank_list <- c(rank_list, format_list[[v_fmt]])
      high_rank <- ifelse(format_list[[v_fmt]] >= min_or_max(rank_list), v_fmt, high_rank)
    }
    message(crayon::yellow(glue::glue("Forcing the metacoder object from the ", crayon::bgWhite(crayon::red({fmt})), " to the ",
                                      crayon::bgWhite(crayon::green({high_rank})),".")))
    mo_clone <- as_MicrobiomeR_format(obj = mo_clone, format = high_rank, ...)
    return(mo_clone)
  } else {
    stop(glue::glue("The metacoder object is not in one of the valid formats: {valid_formats}." ))
  }
}

#' @title As Raw MicrobiomeR Format
#' @description Converts a metacoder object to the raw_format.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @return A Taxmap/metacoder object in the "raw_format".
#' @pretty_print TRUE
#' @details See the [MicrobiomeR_Formats] documentation.
#' @export
#' @family Formatting
#' @rdname as_raw_format
#' @seealso
#'  \code{\link[MicrobiomeR]{is_phyloseq_format}}, \code{\link[MicrobiomeR]{is_raw_format}}, \code{\link[MicrobiomeR]{order_metacoder_data}}
#'  @importFrom crayon silver red greem
as_raw_format <- function(obj) {
  obj <- object_handler(obj = obj)
  mo_clone <- obj$clone()
  if (is_phyloseq_format(mo_clone) == TRUE) {
    mo_clone$data$otu_abundance <- mo_clone$data$otu_table
    mo_clone$data$otu_table <- NULL
    mo_clone$data$otu_annotations <- mo_clone$data$tax_data
    mo_clone$data$tax_data <- NULL
  } else if (is_raw_format(mo_clone)) {
    message(crayon::silver("Converting to the raw format:  The object is already in the raw format."))
  } else {
    stop(crayon::red("Converting to the raw format:  You have to start in the phyloseq format."))
  }
  mo_clone <- order_metacoder_data(obj = mo_clone)
  message(crayon::green("Converted to the raw format."))
  return(mo_clone)
}

#' @title As Basic MicrobiomeR Format
#' @description Converts a metacoder object to the basic_format.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param cols Column names used for \code{\link[metacoder]{calc_taxon_abund}}.  Default: NULL
#' @param out_names Column names of the output used for \code{\link[metacoder]{calc_obs_props}}.  Default: NULL
#' @return A Taxmap/metacoder object in the "basic_format".
#' @pretty_print TRUE
#' @details See the [MicrobiomeR_Formats] documentation.
#' @export
#' @family Formatting
#' @rdname as_basic_format
#' @seealso
#'  \code{\link[MicrobiomeR]{is_phyloseq_format}},  \code{\link[MicrobiomeR]{is_raw_format}},  \code{\link[MicrobiomeR]{is_basic_format}}, \code{\link[MicrobiomeR]{as_raw_format}},  \code{\link[MicrobiomeR]{order_metacoder_data}}
#'
#'  \code{\link[metacoder]{calc_taxon_abund}}, \code{\link[metacoder]{calc_obs_props}}
#' @importFrom metacoder calc_taxon_abund calc_obs_props
#' @importFrom crayon red green silver
as_basic_format <- function(obj, cols = NULL, out_names = NULL) {
  obj <- object_handler(obj = obj)
  mo_clone <- obj$clone()
  # Convert the metacoder object up the heirarchy of formants.
  if (is_phyloseq_format(mo_clone) == TRUE) {
    mo_clone <- as_raw_format(obj = mo_clone)
  }
  # Get metacoder::calc_* parameters
  if (is.null(cols)) {
    cols <- mo_clone$data$sample_data$sample_id
  }
  if (is_raw_format(mo_clone)) {

    suppressMessages({
      suppressWarnings({
        # Create a taxonomy abundance table from the OTU abundance table
        mo_clone$data$taxa_abundance <- metacoder::calc_taxon_abund(obj  = mo_clone,
                                                         data = "otu_abundance",
                                                         cols = cols,
                                                         out_names = out_names)
        # Create an OTU proportions table from the OTU abundance table
        mo_clone$data$otu_proportions <- metacoder::calc_obs_props(obj        = mo_clone,
                                                        data       = "otu_abundance",
                                                        cols       = cols,
                                                        other_cols = TRUE,
                                                        out_names = out_names)
        # Create a taxonomy proportions table from the OTU proportions table
        mo_clone$data$taxa_proportions <- metacoder::calc_taxon_abund(obj  = mo_clone,
                                                           data = "otu_proportions",
                                                           cols = cols,
                                                           out_names = out_names)
        })
      })
  } else if (is_basic_format(mo_clone)) {
    message(crayon::silver("Converting to the basic format:  The object is already in the basic format."))
  } else {
    stop(crayon::red("Converting to the basic format:  You have to start in the phyloseq or raw formats."))
  }
  mo_clone <- order_metacoder_data(obj = mo_clone)
  message(crayon::green("Converted to the basic format."))
  return(mo_clone)
}


#' @title As Analyzed MicrobiomeR Format
#' @description Converts a metacoder object to the analyzed_format.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param cols Column names used for \code{\link[metacoder]{calc_taxon_abund}}.  Default: NULL
#' @param groups Group names used for \code{\link[metacoder]{compare_groups}}.  Default: NULL
#' @param comp_func A Comparison based function used in \code{\link[metacoder]{compare_groups}}.  Default: NULL
#' @param combinations Combinations of treatments used in \code{\link[metacoder]{compare_groups}}.  Default: NULL
#' @param out_names Column names of the output used for \code{\link[metacoder]{calc_obs_props}}.  Default: NULL
#' @return A Taxmap/metacoder object in the "analyzed_format".
#' @pretty_print TRUE
#' @details See the [MicrobiomeR_Formats] documentation.
#' @export
#' @family Formatting
#' @rdname as_analyzed_format
#' @seealso
#'  \code{\link[MicrobiomeR]{is_phyloseq_format}}, \code{\link[MicrobiomeR]{is_raw_format}},  \code{\link[MicrobiomeR]{is_basic_format}},  \code{\link[MicrobiomeR]{is_analyzed_format}},  \code{\link[MicrobiomeR]{as_raw_format}},  \code{\link[MicrobiomeR]{as_basic_format}},  \code{\link[MicrobiomeR]{order_metacoder_data}}
#'
#'  \code{\link[metacoder]{compare_groups}}
#'
#'  \code{\link[taxa]{taxonomy_table}},  \code{\link[taxa]{taxon_ids}}
#' @importFrom metacoder compare_groups
#' @importFrom taxa taxonomy_table taxon_ids
#' @importFrom dplyr rename right_join
#' @importFrom crayon silver red green
as_analyzed_format <- function(obj, cols = NULL, groups = NULL, combinations = NULL, out_names = NULL, comp_func = metacoder_comp_func_1) {
  obj <- object_handler(obj = obj)
  mo_clone <- obj$clone()
  # Convert the metacoder object up the heirarchy of formants.
  if (is_phyloseq_format(mo_clone)) {
    mo_clone <- as_raw_format(obj = mo_clone)
  }
  if (is_raw_format(mo_clone)) {
    mo_clone <- as_basic_format(obj = mo_clone, cols = cols, out_names = out_names)
  }
  # Get metacoder::compare_groups parameters
  if (is.null(cols)) {
    cols <- mo_clone$data$sample_data$sample_id
  }
  if (is.null(groups)) {
    groups <- mo_clone$data$sample_data$TreatmentGroup
  }
  # Continue with conversion to analyzed_format
  if (is_basic_format(mo_clone)) {
    # Compare groups of samples for statistical analysis
    suppressMessages({
      suppressWarnings({
        mo_clone$data$statistical_data <- metacoder::compare_groups(obj = mo_clone,
                                                         data        = "taxa_proportions",
                                                         cols        = cols,
                                                         groups      = groups,
                                                         func        = comp_func,
                                                         other_cols  = TRUE,
                                                         combinations = combinations)
        })
      })
    # Create a table with taxonomy data and stats data for downstream analysis
    tax_table <- obj$taxonomy_table(subset = taxon_ids, add_id_col = TRUE)
    if ("taxon_ids" %in% names(tax_table)) {
      tax_table <- tax_table %>% dplyr::rename(taxon_id = taxon_ids)
    }
    stats_table <- mo_clone$data$statistical_data
    mo_clone$data$stats_tax_data <- dplyr::right_join(x  = tax_table,
                                                      y  = stats_table,
                                                      by = "taxon_id")
  } else if (is_basic_format(mo_clone)) {
    message(crayon::silver("Converting to the analyzed format:  The object is already in the analyzed format."))
  } else {
    stop(crayon::red("Converting to the analyzed format:  You have to start in the phyloseq, raw, or basic formats."))
  }
  # Put data tables in the proper order
  mo_clone <- order_metacoder_data(obj = mo_clone)
  message(crayon::green("Converted to the analyzed format."))
  return(mo_clone)
}


#' @title As MicrobiomeR Format
#' @description Converts a metacoder object to the specified format.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param format The name of the format to convert the object to.
#' @param ... An optional list of parameters to use in the as_*_format function specified
#' by the format parameter.
#' @return A Taxmap/metacoder object in the specified format.
#' @pretty_print TRUE
#' @details See the [MicrobiomeR_Formats] documentation.
#' @export
#' @family Formatting
#' @rdname as_MicrobiomeR_format
#' @seealso
#'  \code{\link[MicrobiomeR]{which_format}},  \code{\link[MicrobiomeR]{as_raw_format}},  \code{\link[MicrobiomeR]{as_basic_format}},  \code{\link[MicrobiomeR]{as_analyzed_format}},  \code{\link[MicrobiomeR]{as_phyloseq_format}},  \code{\link[MicrobiomeR]{object_handler}},  \code{\link[MicrobiomeR]{order_metacoder_data}}
#' @importFrom glue glue
#' @importFrom crayon silver green red
as_MicrobiomeR_format <- function(obj, format, ...) {
  obj <- object_handler(obj = obj)
  mo_clone <- obj$clone()
  current_format <- which_format(mo_clone)
  if (format != current_format) {
    if (format == "raw_format") {
      mo_clone <- as_raw_format(obj = mo_clone)
    } else if (format == "basic_format") {
      mo_clone <- as_basic_format(obj = mo_clone, ...)
    } else if (format == "analyzed_format") {
      mo_clone <- as_analyzed_format(obj = mo_clone, ...)
    } else if (format == "phyloseq_format") {
      mo_clone <- as_phyloseq_format(obj = mo_clone, ...)
    } else {
      stop(crayon::red("The format is not recognized."))
    }
  } else if (format == current_format) {
    message(crayon::silver(glue::glue("Your object is already in the proper format: {format}")))
  } else {
    message(glue::glue("Converted to the ", crayon::green({format}), "."))
  }
  mo_clone <- order_metacoder_data(obj = mo_clone)
  return(mo_clone)
}


#' @title As Phyloseq MicrobiomeR Format
#' @description Converts the metacoder object to the phyloseq_format.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param otu_table The name of the observation table with OTU data.  Default: NULL
#' @param tax_data The name of the observation table with taxonomic annotations.  Default: NULL
#' @param sample_data The name of the observation table with metadata.  Default: NULL
#' @param phy_tree The name of the observation data with the phylogenetic tree.  Default: NULL
#' @return A Taxmap/metacoder object in the phyloseq_format.
#' @details See the [MicrobiomeR_Formats] documentation.
#' @export
#' @family Formatting
#' @rdname as_phyloseq_format
#' @seealso
#'  \code{\link[MicrobiomeR]{object_handler}},\code{\link[MicrobiomeR]{order_metacoder_data}}
as_phyloseq_format <- function(obj, otu_table="otu_abundance", tax_data="otu_annotations", sample_data="sample_data", phy_tree="phy_tree") {
  obj <- object_handler(obj = obj)
  mo_clone <- obj$clone()
  if (!is.null(otu_table)) {
    mo_clone$data$otu_table <- mo_clone$data[otu_table]
    mo_clone$data[otu_table] <- NULL
  }
  if (!is.null(tax_data)) {
    mo_clone$data$tax_data <- mo_clone$data[tax_data]
    mo_clone$data[tax_data] <- NULL
  }
  if (!is.null(sample_data)) {
    mo_clone$data$sample_data <- mo_clone$data[sample_data]
    mo_clone$data[sample_data] <- NULL
  }
  if (!is.null(phy_tree)) {
    mo_clone$data$phy_tree <- mo_clone$data[phy_tree]
    mo_clone$data[phy_tree] <- NULL
  }
  mo_clone <- order_metacoder_data(obj = mo_clone)
  return(mo_clone)
}

#' @title As Custom MicrobiomeR Format
#' @description A function for formatting metacoder objects in the MicrobiomeR format.  This function
#' attempts to give more customization than the as_*_format functions.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param format The name of the format to convert the object to.
#' @param change_name_list A list with names of the tables in the observation data, that have matching
#' values that are used to change the names of the table.   Default: NULL
#' @param ... An optional list of parameters to use in the as_*_format function specified
#' @return A metacoder object that we have tried to format with all of our heart.
#' @pretty_print TRUE
#' @details This function is meant to be more helpful for customizing the metacoder object.
#' @export
#' @family Formatting
#' @rdname as_custom_format
#' @seealso
#'  \code{\link[MicrobiomeR]{object_handler}},  \code{\link[MicrobiomeR]{which_format}},  \code{\link[MicrobiomeR]{order_metacoder_data}},  \code{\link[MicrobiomeR]{as_MicrobiomeR_format}}
#' @importFrom glue glue
#' @importFrom crayon red
as_custom_format <- function(obj, format, change_name_list = NULL, ...) {

  # Metacoder Objects
  obj <- object_handler(obj = obj)
  mo_clone <- obj$clone()
  fmt <- which_format(obj = mo_clone)
  # Get table lists and vectors
  format_level_list <- pkg.private$format_level_list

  # Create vector of current table names
  obs_tables <- names(mo_clone$data)

  # Logic for getting to the right format
  if (fmt == format) {
    mo_clone <- order_metacoder_data(obj = mo_clone)
    return(mo_clone)
  } else if (!is.null(change_name_list)) {   # Create a list of tables to change if necessary
    # Create a list of key/values used to change names
    changed_tables <- change_name_list[names(change_name_list) %in% obs_tables]
    # Create a list of key/values used to create new tables
    bad_table_names <- change_name_list[!names(change_name_list) %in% obs_tables]

    # Throw errors for bad table names
    if (length(bad_table_names != 0)) {
      if (length(bad_table_names) == length(change_name_list)) {
        stop(glue::glue(crayon::red("None of the parameters that you've given are in your observation data:
             {bad_table_names}")))
      } else {
        stop(glue::glue(crayon::red("You have given some bad table names that aren't in you metacoder object:
                       {bad_table_names}")))
      }
    }
    # Change the table names
    for (current_table in names(changed_tables)) {
      # Create a new table from the data in the current table
      mo_clone$data[changed_tables[[current_table]]] <- mo_clone$data[[current_table]]
      # Remove the current table from the observation list
      mo_clone$data[current_table] <- NULL
    }
    # See if the changed names helped with the format
    changed_obs_tables <- names(mo_clone$data)
    fmt = which_format(obj = mo_clone)
    if (fmt == format) {
      mo_clone <- order_metacoder_data(obj = mo_clone)
      return(mo_clone)
    }
  }
  # If the current format is at a lower level than the expected format continue
  if (format_level_list[[fmt]] < format_level_list[[format]]) {
    # If the level is not negative continue
    if (sign(format_level_list[[fmt]]) == 1) {
      mo_clone <- as_MicrobiomeR_format(obj = mo_clone, format = format, ...)
    } else { # Throw an error if the level is negative (unknown or mixed format)
      warning(glue::glue(crayon::red("Here is a list of your observation data:
                           {changed_obs_tables}")))
      stop(crayon::red("Your data is in an unknown or mixed format."))
    }
  }
  mo_clone <- order_metacoder_data(obj = mo_clone)
  return(mo_clone)
}



