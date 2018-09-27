#' @title Get Alpha Diversity
#'
#' @description This function generates various alpha diversity measures include Shannon, Fisher, Coverage, Gini Simpson, and Inverse Simpson.
#'
#' @param phyloseq_object A preprocessed phyloseq object.
#'
#' @return Returns a list object of alpha diversity measurements.
#' @export
get_alpha_diversity <- function(phyloseq_object) {
  # Get all of the diversities.
  divs <- microbiome::diversities(phyloseq_object, index = "all")
  phyloseq_object.meta <- microbiome::meta(phyloseq_object)
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
