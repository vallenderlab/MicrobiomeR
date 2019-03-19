#' @title Alpha Diversity Measures
#' @description This function generates various alpha diversity measures include Shannon, Fisher, Coverage, GiniSimpson, and InverseSimpson.
#' @param obj An object to be converted to a Taxmap object with \code{\link[MicrobiomeR]{create_taxmap}}.
#' @param group The "TreatmentGroup" or similar grouping from your metadata to denote sample groups, Default: 'TreatmentGroup'
#' @return Returns a list of alpha diversity measures with metadata.
#' @pretty_print TRUE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(MicrobiomeR)
#'   data <- analyzed_silva
#'   measures <- alpha_diversity_measures(data)
#'   measures$Shannon
#' }
#' }
#' @export
#' @rdname alpha_diversity_measures
#' @seealso
#'  View \code{\link[microbiome]{alpha}},\code{\link[microbiome]{meta}} to learn more about generating diversity measures with the microbiome package.
#' @importFrom microbiome alpha meta
#' @importFrom metacoder as_phyloseq
#' @importFrom utils combn
alpha_diversity_measures <- function(obj, group = "TreatmentGroup") {
  metacoder_object <- validate_MicrobiomeR_format(
    obj = create_taxmap(obj),
    valid_formats = c("analyzed_format")
  )

  # Convert Taxmap object to a phyloseq object.
  phyloseq_object <- metacoder::as_phyloseq(metacoder_object, otu_table = "otu_abundance", phy_tree = "phy_tree")

  # Get all of the diversities.
  divs <- microbiome::alpha(phyloseq_object, index = "all")
  phyloseq_object.meta <- microbiome::meta(phyloseq_object)
  phyloseq_object.meta$Shannon <- divs$diversity_shannon
  phyloseq_object.meta$InverseSimpson <- divs$diversity_inverse_simpson
  phyloseq_object.meta$GiniSimpson <- divs$diversity_gini_simpson
  phyloseq_object.meta$Fisher <- divs$diversity_fisher
  phyloseq_object.meta$Coverage <- divs$diversity_coverage

  # create a list of pairwise comaprisons
  groups <- levels(as.factor(phyloseq_object.meta[[group]])) # get the variables
  num_groups <- length(groups)

  # make a pairwise list that we want to compare.
  group.pairs <- utils::combn(seq_along(groups), num_groups, simplify = FALSE, FUN = function(i) groups[i])

  phyloseq_object.meta$group.pairs <- group.pairs

  return(phyloseq_object.meta)
}

#' @title Alpha Diversity Plot
#' @description Plot the alpha diversity using a violin plot. `alpha_diversity_plots` generates plots for all alpha diversity measures.
#' @param obj An object to be converted to a Taxmap object with \code{\link[MicrobiomeR]{create_taxmap}}.
#' @param measure Select an alpha diversity measure such as Shannon, Fisher, Coverage, GiniSimpson, and InverseSimpson, Default: 'Shannon'
#' @param group The "TreatmentGroup" or similar grouping or column from your metadata to denote sample groups, Default: 'TreatmentGroup'
#' @param select_otu_table Choose an otu table to analyze, Default: 'otu_proportions'
#' @param title The title of the plot, Default: NULL
#' @return Returns an alpha diversity plot.
#' @details Alpha diversity helps to determine the species richness (the number of different species in a sample) or evenness (similar abundance level).
#' We prefer to use `shannon` as it is better for data generated using the QIIME pipeline.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(MicrobiomeR)
#'   data <- analyzed_silva
#'   plot <- alpha_diversity_plot(obj = data,
#'                                measure = "Shannon",
#'                                select_otu_table = "otu_proportions")
#'   plot
#' }
#' }
#' @export
#' @rdname alpha_diversity_plot
#' @seealso \code{\link{alpha_diversity_measures}}, \code{\link[vegan]{diversity}}, \code{\link[ggpubr]{ggviolin}}
#' @family Visualizations
#' @importFrom ggplot2 xlab ylab
#' @importFrom ggpubr stat_compare_means ggviolin
#' @importFrom ggthemes theme_pander
#' @importFrom utils combn
#' @importFrom vegan diversity
alpha_diversity_plot <- function(obj, measure = "Shannon", group = "TreatmentGroup", title = NULL) {
  # Validate data format
  metacoder_object <- validate_MicrobiomeR_format(
    obj = create_taxmap(obj),
    valid_formats = c("analyzed_format")
  )
  measures <- alpha_diversity_measures(obj = metacoder_object, group = group)
  metacoder_object$data$sample_data[[measure]] <- measures[[measure]]

  if (typeof(metacoder_object$data$sample_data$TreatmentGroup) == "character") {
    metacoder_object$data$sample_data[[group]] <- as.factor(metacoder_object$data$sample_data$TreatmentGroup)
  } else {
    metacoder_object$data$sample_data[[group]] <- factor(metacoder_object$data$sample_data[[group]], levels = list(unique(metacoder_object$data$sample_data[[group]])))
  }
  groups <- levels(metacoder_object$data$sample_data[[group]]) # get the variables
  num_groups <- length(groups)

  # make a pairwise list that we want to compare.
  metacoder_object$data$sample_data$group.pairs <- utils::combn(seq_along(groups), num_groups, simplify = FALSE, FUN = function(i) groups[i])

  if (num_groups == 2) {
    palette <- c("#3288bd", "#d53e4f")
  } else {
    # Create a palette when there are more than 2 groups
    palette <- c("#3288bd", "#d53e4f", "#62954C", "#C59144")
  }

  plot <- ggpubr::ggviolin(metacoder_object$data$sample_data,
    x = group, y = measure,
    color = "black",
    add = "boxplot",
    fill = group,
    palette = palette,
    legend.title = title
  ) + ggplot2::xlab(title) + ggplot2::ylab(toupper(measure)) +
    ggthemes::theme_pander() + ggpubr::stat_compare_means(comparisons = metacoder_object$data$sample_data$group.pairs, label = "p.signif", label.y = 7) +
    ggpubr::stat_compare_means(label.y = 8)

  return(plot)
}

#' @param measures A list of alpha diversity measures such as Shannon, Fisher, Coverage, GiniSimpson, and InverseSimpson, Default: 'c("Shannon", "GiniSimpson", "InverseSimpson")'
#' @return Returns a melted dataframe.
#' @family Visualizations
#' @rdname alpha_diversity_plot
alpha_diversity_plots <- function(obj, measures = c("Shannon", "GiniSimpson", "InverseSimpson"), group = "TreatmentGroup") {
  if (is.null(measures)) {
    measures <- c("Shannon", "GiniSimpson", "InverseSimpson")
  } else if (length(measures) < 2) {
    stop("Use the alpha_diversity_plot function for generating a plot for 1 alpha diversity index.")
  }
  alpha_div_plots <- list()
  for (m in measures) {
    alpha_div_plots[[m]] <- alpha_diversity_plot(obj, measure = m, group = group, title = m)
  }

  return(alpha_div_plots)
}

#' @title Save Alpha Diversity Plots
#' @description This function saves alpha diversity plots stored in a list object to an output folder.
#' @param alpha_div_plots A named list of alpha diversity plots.
#' @param format The format of the output image.  Default: 'tiff'
#' @param start_path The starting path of the output directory.  Default: 'output'
#' @param ... An optional list of parameters to use in the output_dir function.
#' @return An output directory that contains alpha diversity plots.
#' @pretty_print TRUE
#' @details This function creates an appropriate output directory, where it saves publication ready
#' plots.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # This example uses data that are no longer available in the MicrobiomeR package,
#'   # however, they can be easily generated with \code{\link{MicrobiomeR}{as_analyzed_format}}.
#'   library(MicrobiomeR)
#'   analyzed_silva <- as_MicrobiomeR_format(MicrobiomeR::raw_silva_2, "analyzed_format")
#'   alpha_div_plots <- alpha_diversity_plots(analyzed_silva)
#'   # Save to \emph{./output/alpha_diversity} folder.
#'   save_alpha_diversity_plot(alpha_div_plots)
#' }
#' }
#' @export
#' @family Visualizations
#' @rdname save_alpha_diversity_plots
#' @seealso
#'  \code{\link[MicrobiomeR]{output_dir}}
#'
#'  \code{\link[ggplot2]{ggsave}}
#' @importFrom ggplot2 ggsave
#' @importFrom crayon yellow green
#' @importFrom glue glue
save_alpha_diversity_plots <- function(alpha_div_plots, format = "tiff", start_path = "output", ...) {
  # Create the relative path to the alpha diversity plots.  By default the path will be <pwd>/output/<experiment>/alpha_diversity/<format(Sys.time(), "%Y-%m-%d_%s")>
  # With the parameters set the full path will be <pwd>/output/<experiment>/alpha_diversity/<extra_path>.
  full_path <- output_dir(start_path = start_path, plot_type = "alpha_diversity", ...)
  message(glue::glue(crayon::yellow("Saving Alpha Diversity plots to the following directory: \n", "\r\t{full_path}")))
  # Iterate the plot list and save them in the proper directory
  for (measure in names(alpha_div_plots)) {
    if (measure != "metacoder_object") {
      message(crayon::green("Saving the {measure} alpha aiversity plot."))
      ggplot2::ggsave(filename = sprintf("%s_alpha_diversity.%s", measure, format), plot = alpha_div_plots[[measure]], device = format, path = full_path, dpi = 500)
    }
  }
}
