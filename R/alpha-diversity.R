#' @title Get Alpha Diversity
#' @description This function generates various alpha diversity measures include Shannon, Fisher, Coverage, Gini Simpson, and Inverse Simpson.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
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
  metacoder_object <- validate_MicrobiomeR_format(
    obj = object_handler(obj),
    valid_formats = c("analyzed_format")
  )

  # Convert metacoder object to a phyloseq object.
  phyloseq_object <- as_phyloseq_format(metacoder_object)

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

alpha_diversity_plot <- function(obj, measure = "Shannon", select_otu_table = "otu_proportions", save = FALSE) {
  # Validated data format
  metacoder_object <- validate_MicrobiomeR_format(
    obj = object_handler(obj),
    valid_formats = c("analyzed_format")
  )
  metacoder_object$data$sample_data[[measure]] <- vegan::diversity(metacoder_object$data[[select_otu_table]][, metacoder_object$data$sample_data$X.SampleID],
    MARGIN = 2, index = measure
  )

  metacoder_object$data$sample_data$TreatmentGroup <- factor(metacoder_object$data$sample_data$TreatmentGroup, levels = c("Control", "Stressed"))
  treatment_groups <- levels(metacoder_object$data$sample_data$TreatmentGroup) # get the variables

  # make a pairwise list that we want to compare.
  metacoder_obj$data$sample_data$treatment_group.pairs <- combn(seq_along(treatment_groups), 2, simplify = FALSE, FUN = function(i) treatment_groups[i])

  plot <- ggviolin(metacoder_object$data$sample_data,
    x = "TreatmentGroup", y = measure,
    color = "black",
    add = "boxplot",
    fill = "TreatmentGroup",
    palette = c("#3288bd", "#d53e4f"),
    legend.title = "Treatment Group"
  ) + xlab("Treatment Group") + ylab(toupper(measure)) +
    theme_pander() + stat_compare_means(comparisons = metacoder_object$data$sample_data$treatment_group.pairs, label = "p.signif", label.y = 7) +
    stat_compare_means(label.y = 8)

  if (save == TRUE) {
    ggsave(filename = "alpha_diversity.tiff")
    return(plot)
  } else {
    return(plot)
  }
}
