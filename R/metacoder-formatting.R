# TODO: Continue to update the pure documentation.  basic/analyzed formats; dataManip/Visualizations

#' The MicrobiomeR Formats
#'
#' The formats are listed in order from lowest to highest level.
#' To get a taxmap object into the higher formats they must be in one of the lower formats.
#' It is also inherently important to understand the [MicrobiomeR_Workflow], which describes some of the important
#' dependencies used to create the MicrobiomeR formats listed below:
#' \describe{
#'   \item{phyloseq_format}{Level 0.  The "phyloseq_format" is the format given to a \emph{taxmap} object that has just
#'   been converted from a \emph{phyloseq} object using \emph{metacoder} via \code{\link[MicrobiomeR]{get_phyloseq_obj}}.  A taxmap object identified by MicrobiomeR
#'   to be in this format will contain "otu_table", "tax_data", "sample_data" and "phy_tree" data in the observation
#'   tables (\strong{metacoder_obj$data}), which represents the matching \emph{phyloseq::otu_table}, \emph{phyloseq::tax_table},
#'   \emph{phyloseq::sample_data}, and \emph{phyloseq::phy_tree} data in the original phyloseq object.}
#'   \item{raw_format}{Level 1.  The "raw_format" is the format given to a \emph{taxmap} object that has been
#'   processed with \code{\link[MicrobiomeR]{as_raw_format}}.  It contains the "otu_abundance" and "otu_annotations" data in the observation
#'   tables, which are just name conversions of the \emph{otu_table} and \emph{tax_data} table from the "phyloseq_format".}
#'   \item{basic_format}{Level 2.  The "basic_format" is the format given to a \emph{taxmap} object that has been
#'   processed with \code{\link[MicrobiomeR]{as_basic_format}}.  It contains the tables available in the "raw_format" along with
#'   "taxa_abundance", "otu_proportions", and "taxa_proportions" tables.}
#'   \item{analyzed_format}{Level 3.  The "analyzed_format" is the format given to a \emph{taxmap} object that has been
#'   processed with \code{\link[MicrobiomeR]{as_analyzed_format}}.  It contains the tables available in the "basic_format" along
#'   with "statistical_data" and "stats_tax_data" tables.}
#' }
#'
#' @name MicrobiomeR_Formats
NULL

#' The MicrobiomeR Workflow
#'
#' In order to understand the MicrobiomeR packages workflow it is necessary to understand a few other
#' packages.  Primarily \strong{phyloseq}, \strong{taxa}, and \strong{metacoder}.
#'
#' @section Data Import:
#' Data is imported with the \code{\link[MicrobiomeR]{get_phyloseq_obj}} function or something proprietary that
#' creates a \code{\link[phyloseq:phyloseq-class]{phyloseq object}} with an \code{\link[phyloseq:otu_table-class]{phyloseq::otu_table}},
#' a \code{\link[phyloseq:taxonomyTable-class]{phyloseq::tax_table}}, a \code{\link[phyloseq:sample_data-class]{phyloseq::sample_data}},
#' and a \code{\link[phyloseq:phylo-class]{phyloseq::phy_tree}}.  The data is preferrably not filtered
#' with phyloseq at this stage.
#' @section Data Manipulation:
#' Data is manipulated into a \code{\link[taxa:taxmap]{taxa::taxmap}} object using the \code{\link[metacoder:parse_phyloseq]{metacoder::parse_phyloseq}}
#' function.  Custom [MicrobiomeR_Formats] are used to validate the data for consistency.  MicrobiomeR and metacoder filtering
#' functions can be used to preprocess the data before visualization.
#' @section Visualizations:
#' MicrobiomeR contains useful functions for creating visualizations such as heat_trees, correlations plots, stacked bar charts,
#' violin plots, and PCA plots.  These functions are also coupled with save functions that organize saved output by experiment.
#'
#' @name MicrobiomeR_Workflow
NULL

#' @title Which MicrobiomeR Format
#' @description A function for looking at a metacoder object and returning the identified MicrobiomeR format.
#' @param obj A Taxmap/metacoder object.
#' @return If the format is verified it returns a character string denoting the identified format.
#' @pretty_print TRUE
#' @details This function is used to get basic infomration about the format of the taxmap object
#' that is supplied by the user.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Formatting and Validation
#' @rdname which_format
#' @seealso
#'   \code{\link[MicrobiomeR]{is_raw_format}},  \code{\link[MicrobiomeR]{is_basic_format}},  \code{\link[MicrobiomeR]{is_analyzed_format}},  \code{\link[MicrobiomeR]{is_phyloseq_format}}
#'   \code{\link[MicrobiomeR]{as_raw_format}},  \code{\link[MicrobiomeR]{as_basic_format}},  \code{\link[MicrobiomeR]{as_analyzed_format}},  \code{\link[MicrobiomeR]{as_phyloseq_format}}
#'   \code{\link[MicrobiomeR]{as_MicrobiomeR_format}},  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}}
which_format <- function(obj) {
  mo_clone <- obj$clone()
  # Table names
  table_names <- names(mo_clone$data)
  # Format names
  raw_names <- c("otu_abundance", "otu_annotations")
  basic_names <- c(raw_names, "taxa_abundance", "otu_proportions", "taxa_proportions")
  analyzed_names <- c(basic_names, "statistical_data", "stats_tax_data")
  phyloseq_names <- c("otu_table", "tax_data", "sample_data", "phy_tree")
  # Flags
  raw_flag <- all(raw_names %in% table_names)
  basic_flag <- all(basic_names %in% table_names)
  analyzed_flag <- all(analyzed_names %in% table_names)
  phyloseq_flag <- all(phyloseq_names %in% table_names)
  if (phyloseq_flag) {
    other_flag <- any(c(raw_flag, basic_flag, analyzed_flag))
    if (!other_flag) {
      warning("Your object is in the phyloseq format!")
      warning("Please format your metacoder object to continue analysis.")
      return("phyloseq_format")
    } else {
      warning(sprintf("The table names in the metacoder object are: %s", paste(table_names, collapse = ", ")))
      stop("You have a mix between phyloseq format and other format.")
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
      warning(sprintf("The table names in the metacoder object are: %s", paste(table_names, collapse = ", ")))
      stop("The object is not in a recognized format.")
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
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Formatting and Validation
#' @rdname is_raw_format
#' @seealso
#'  \code{\link[MicrobiomeR]{which_format}}
is_raw_format <- function(obj) {
  fmt <- MicrobiomeR::which_format(obj)
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
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Formatting and Validation
#' @rdname is_basic_format
#' @seealso
#'  \code{\link[MicrobiomeR]{which_format}}
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
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Formatting and Validation
#' @rdname is_analyzed_format
#' @seealso
#'  \code{\link[MicrobiomeR]{which_format}}
is_analyzed_format <- function(obj) {
  fmt <- which_format(obj)
  if (fmt == "analyzed_format"){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# Returns TRUE if the MicrobiomeR format is the phyloseq_format
#' @title Is Phyloseq MicrobiomeR Format
#' @description This function returns a logical based on weather or not the object is in the phyloseq_format.
#' @param obj A Taxmap/metacoder object.
#' @return A logical (TRUE/FALSE).
#' @pretty_print TRUE
#' @details The "phyloseq_format" is Level 0. in the [MicrobiomeR_Formats] hierarchy.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Formatting and Validation
#' @rdname is_phyloseq_format
#' @seealso
#'  \code{\link[MicrobiomeR]{which_format}}
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
#' @param obj A Taxmap/metacoder object.
#' @return A taxmap object with observation data in the proper order for downstream analysis.
#' @pretty_print TRUE
#' @details Changes the order of the observation tables in \strong{metacoder_object$data} to
#' otu_abundance, otu_annotations, otu_proportions, sample_data, phy_tree, taxa_abundance,
#' taxa_proportions, statistical_data, and stats_tax_data respectively.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Formatting and Validation
#' @rdname order_metacoder_data
order_metacoder_data <- function(obj) {
  mo_clone <- obj$clone()
  # Create list of expected table names
  expected_names <- list(otu_abundance    = "otu_abundance",
                         otu_annotations  = "otu_annotations",
                         otu_proportions  = "otu_proportions",
                         sample_data      = "sample_data",
                         phy_tree         = "phy_tree",
                         taxa_abundance   = "taxa_abundance",
                         taxa_proportions = "taxa_proportions",
                         statistical_data = "statistical_data",
                         stats_tax_data   = "stats_tax_data")

  # Create vector that contains expected tables names already in the metacoder object
  table_order <- names(expected_names[names(expected_names) %in% names(mo_clone$data)])
  # Get the table names that aren't expected but in the metacoder object
  other_names <- names(mo_clone$data)[!names(mo_clone$data) %in% c(expected_names)]
  table_order <- c(table_order, other_names)
  mo_clone$data <- mo_clone$data[table_order]
  return(mo_clone)
}
