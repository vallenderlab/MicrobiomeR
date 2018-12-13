#' @title Get Alpha Diversity
#' @description This function generates various alpha diversity measures include Shannon, Fisher, Coverage, Gini Simpson, and Inverse Simpson.
#' @param obj A phyloseq or metacoder object.
#' @return Returns a list an alpha diversity object.
#' @pretty_print TRUE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname get_alpha_diversity
#' @seealso
#'  \code{\link[microbiome]{diversities}},\code{\link[microbiome]{meta}}
#' @importFrom microbiome diversities meta
get_alpha_diversity <- function(obj) {
  # Get all of the diversities.
  divs <- microbiome::diversities(obj, index = "all")
  phyloseq_object.meta <- microbiome::meta(obj)
  phyloseq_object.meta$Shannon <- divs$shannon
  phyloseq_object.meta$InverseSimpson <- divs$inverse_simpson
  phyloseq_object.meta$GiniSimpson <- divs$gini_simpson
  phyloseq_object.meta$Fisher <- divs$fisher
  phyloseq_object.meta$Coverage <- divs$coverage

  # create a list of pairwise comaprisons
  groups <- levels(phyloseq_object.meta$TreatmentGroup) # get the variables

  # make a pairwise list that we want to compare.
  group.pairs <- combn(seq_along(groups), 2, simplify = FALSE, FUN = function(i) groups[i])

  phyloseq_object.meta$group.pairs <- group.pairs

  return(phyloseq_object.meta)
}
