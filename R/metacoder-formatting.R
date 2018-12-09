#' @title Which MicrobiomeR Format
#' @description A function for looking at a metacoder object and returning the identified MicrobiomeR format.
#' @param obj A Taxmap/metacoder object
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
