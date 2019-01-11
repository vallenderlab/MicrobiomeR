#' @title Get Alpha Diversity
#' @description This function generates various alpha diversity measures include Shannon, Fisher, Coverage, Gini Simpson, and Inverse Simpson.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param group The "TreatmentGroup" or similar grouping from your metadata to denote sample groups, Default: 'TreatmentGroup'
#' @return Returns a list of alpha diversity measures with metadata.
#' @pretty_print TRUE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(MicrobiomeR)
#'   data <- analyzed_silva
#'   measures <- get_alpha_diversity_measures(data)
#'   measures$Shannon
#' }
#' }
#' @export
#' @rdname get_alpha_diversity_measures
#' @seealso
#'  View \code{\link[microbiome]{diversities}},\code{\link[microbiome]{meta}} to learn more about generating diversity measures with the microbiome package.
#' @importFrom microbiome diversities meta
#' @importFrom metacoder as_phyloseq
#' @importFrom utils combn
get_alpha_diversity_measures <- function(obj, group = "TreatmentGroup") {
  metacoder_object <- validate_MicrobiomeR_format(
    obj = object_handler(obj),
    valid_formats = c("analyzed_format")
  )

  # Convert metacoder object to a phyloseq object.
  phyloseq_object <- metacoder::as_phyloseq(metacoder_object, otu_table = "otu_abundance", phy_tree = "phy_tree")

  # Get all of the diversities.
  divs <- microbiome::diversities(phyloseq_object, index = "all")
  phyloseq_object.meta <- microbiome::meta(phyloseq_object)
  phyloseq_object.meta$Shannon <- divs$shannon
  phyloseq_object.meta$InverseSimpson <- divs$inverse_simpson
  phyloseq_object.meta$GiniSimpson <- divs$gini_simpson
  phyloseq_object.meta$Fisher <- divs$fisher
  phyloseq_object.meta$Coverage <- divs$coverage

  # create a list of pairwise comaprisons
  groups <- levels(as.factor(phyloseq_object.meta[[group]])) # get the variables

  # make a pairwise list that we want to compare.
  group.pairs <- utils::combn(seq_along(groups), 2, simplify = FALSE, FUN = function(i) groups[i])

  phyloseq_object.meta$group.pairs <- group.pairs

  return(phyloseq_object.meta)
}

#' @title Alpha Diversity Plot
#' @description Plot the alpha diversity using a violin plot.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param measure Select an alpha diversity measure such as shannon, gini simpson, and inverse simpson, Default: 'shannon'
#' @param select_otu_table Choose an otu table to analyze, Default: 'otu_proportions'
#' @param save Save the plot, Default: FALSE
#' @return Returns an alpha diversity plot.
#' @details Alpha diversity helps to determine the species richness (the number of different species in a sample) or evenness (similar abundance level).
#' We prefer to use `shannon` as it is better for data generated using the QUIIME pipeline.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(MicrobiomeR)
#'   data <- analyzed_silva
#'   plot <- alpha_diversity_plot(obj = data, measure = "shannon", select_otu_table = "otu_proportions", save = FALSE)
#'   plot
#' }
#' }
#' @export
#' @rdname alpha_diversity_plot
#' @seealso \code{\link{get_alpha_diversity_measures}}, \code{\link[vegan]{diversity}}, \code{\link[ggpubr]{ggviolin}}
#' @family Visualizations
#' @importFrom ggplot2 xlab ylab
#' @importFrom ggpubr stat_compare_means ggviolin
#' @importFrom ggthemes theme_pander
#' @importFrom utils combn
alpha_diversity_plot <- function(obj, measure = "shannon", select_otu_table = "otu_proportions") {
  # Validate data format
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
  metacoder_object$data$sample_data$treatment_group.pairs <- utils::combn(seq_along(treatment_groups), 2, simplify = FALSE, FUN = function(i) treatment_groups[i])

  plot <- ggpubr::ggviolin(metacoder_object$data$sample_data,
    x = "TreatmentGroup", y = measure,
    color = "black",
    add = "boxplot",
    fill = "TreatmentGroup",
    palette = c("#3288bd", "#d53e4f"),
    legend.title = "Treatment Group"
  ) + ggplot2::xlab("Treatment Group") + ggplot2::ylab(toupper(measure)) +
    ggthemes::theme_pander() + ggpubr::stat_compare_means(comparisons = metacoder_object$data$sample_data$treatment_group.pairs, label = "p.signif", label.y = 7) +
    ggpubr::stat_compare_means(label.y = 8)

  return(plot)
}
