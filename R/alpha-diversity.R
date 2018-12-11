#' @title Get Alpha Diversity
#'
#' @description This function generates various alpha diversity measures include Shannon, Fisher, Coverage, Gini Simpson, and Inverse Simpson.
#'
#' @param phyloseq_object A phyloseq or metacoder object.
#'
#' @return Returns a list object of alpha diversity measurements.
#' @importFrom microbiome diversities meta
#' @importFrom utils combn
#' @export
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
