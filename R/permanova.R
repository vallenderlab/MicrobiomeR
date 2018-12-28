#' @title PERMANOVA significance test for group-level differences
#' @description Evaluate whether the group(s) has a significant effect on overall gut microbiota composition.
#' @param obj Import an unprocessed phyloseq object.
#' @param distance_method Use the desired statistical method. The default method is "bray".
#' @return Returns a list object of the data.
#' @export
permanova <- function(obj, distance_method = "bray") {
  # Validate data format
  metacoder_object <- validate_MicrobiomeR_format(
    obj = object_handler(obj),
    valid_formats = c("analyzed_format")
  )
  # TODO: Convert this to metacoder strictly.
  phyloseq <- list()
  otu <- microbiome::abundances(phyloseq_object)
  meta <- microbiome::meta(phyloseq_object)
  if (distance_method == "wunifrac" | distance_method == "unifrac") {
    stop("Method is not integrated yet.")
  } else {
    phyloseq[["permanova"]] <- vegan::adonis(t(otu) ~ TreatmentGroup, data = meta, permutations = 99, method = distance_method)
  }
  phyloseq[["phyloseq"]] <- phyloseq_object

  # Checking the homogeneity condition
  # Note the assumption of similar multivariate spread among the groups
  # ie. analogous to variance homogeneity
  # Here the groups have signif. different spreads and
  # permanova result may be potentially explained by that.
  dist <- vegan::vegdist(t(otu))
  phyloseq[["anova"]] <- anova(vegan::betadisper(dist, meta$TreatmentGroup))

  # Investigate the top factors
  # Show coefficients for the top taxa separating the groups
  phyloseq[["coef"]] <- coefficients(phyloseq$permanova)["TreatmentGroup1", ]
  phyloseq[["top.coef"]] <- phyloseq$coef[rev(order(abs(phyloseq$coef)))[1:50]]
  return(phyloseq)
}

#' @title Top Coefficients Barplot
#' @description FUNCTION_DESCRIPTION
#' @param top_coefficients PARAM_DESCRIPTION
#' @return Returns a barplot of the top coefficients
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname topcoef_barplot
topcoef_barplot <- function(top_coefficients) {
  # Set graphical parameters
  par(mar = c(3, 16, 2, 2))

  # Plot the top coefficients
  plot <- barplot(sort(top_coefficients),
    horiz = T, las = 1, main = "Top taxa",
    col = ifelse(sort(top_coefficients) >= 0, "#3288bd", "#d53e4f")
  )
  return(plot)
}
