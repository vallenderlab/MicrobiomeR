# TODO: Continue to update the pure documentation.  basic/analyzed formats; dataManip/Visualizations

#' The MicrobiomeR Formats
#'
#' The formats are listed in order from lowest to highest level.
#' To get a taxmap object into the higher formats they must be in one of the lower formats.
#' It is also inherently important to understand the [MicrobiomeR_Workflow], which describes some of the important
#' dependencies used to create the MicrobiomeR formats listed below:
#' \describe{
#'   \item{phyloseq_format}{Level 0.  The "phyloseq_format" is the format given to a \emph{taxmap} object that has just
#'   been converted from a \emph{phyloseq} object using \emph{metacoder} via \code{\link[MicrobiomeR]{create_phyloseq}}.  A taxmap object identified by MicrobiomeR
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
#' @family Formatting
#' @name MicrobiomeR_Formats
NULL

#' The MicrobiomeR Workflow
#'
#' In order to understand the MicrobiomeR packages workflow it is necessary to understand a few other
#' packages.  Primarily \strong{phyloseq}, \strong{taxa}, and \strong{metacoder}.
#'
#' @section Data Import:
#' Data is imported with the \code{\link[MicrobiomeR]{create_phyloseq}} function or something proprietary that
#' creates a \code{\link[phyloseq:phyloseq-class]{phyloseq object}} with an \code{\link[phyloseq:otu_table-class]{phyloseq::otu_table}},
#' a \code{\link[phyloseq:taxonomyTable-class]{phyloseq::tax_table}}, a \code{\link[phyloseq:sample_data-class]{phyloseq::sample_data}},
#' and a \code{\link[phyloseq:phylo-class]{phyloseq::phy_tree}}.  The data is preferrably not filtered
#' with phyloseq at this stage.
#'
#' @section Data Manipulation:
#' Data is manipulated into a \code{\link[taxa:taxmap]{taxa::taxmap}} object using the \code{\link[metacoder:parse_phyloseq]{metacoder::parse_phyloseq}}
#' function.  Custom [MicrobiomeR_Formats] are used to validate the data for consistency.  MicrobiomeR and metacoder filtering
#' functions can be used to preprocess the data before visualization.
#'
#' @section Visualizations:
#' MicrobiomeR contains useful functions for creating visualizations such as heat_trees, correlations plots, stacked bar charts,
#' violin plots, and PCA plots.  These functions are also coupled with save functions that organize saved output by experiment.
#'
#' @name MicrobiomeR_Workflow
NULL
