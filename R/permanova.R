#' @title PERMANOVA
#' @description Evaluate whether the group(s) has a significant effect on overall gut microbiota composition.
#' @param obj An object to be converted to a taxmap object with \code{\link[MicrobiomeR]{create_taxmap}}.
#' @param distance_method  Use a desired distance method, Default: 'bray'
#' @param group The group or column in the metadata to test upon, Default: 'TreatmentGroup'
#' @return Returns a list which includes permanova, anova, coefficients, and top coefficients.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(MicrobiomeR)
#'   data <- analyzed_silva
#'   p <- permanova(obj = data, group = "TreatmentGroup")
#'   p$permanova
#' }
#' }
#' @export
#' @rdname permanova
#' @seealso View \code{\link{top_coefficients_barplot}} to plot the top_coefficients returned from this function.
#' See \code{\link[vegan]{adonis}} and \code{\link[vegan]{betadisper}} to understand more about how the permanova data was generated
#' @importFrom dplyr select
#' @importFrom microbiome abundances meta
#' @importFrom vegan adonis vegdist betadisper
#' @importFrom phyloseq distance
#' @importFrom dplyr expr sym enquo
permanova <- function(obj, distance_method = "bray", group = "TreatmentGroup") {
  # Validate data format
  metacoder_object <- validate_MicrobiomeR_format(
    obj = create_taxmap(obj),
    valid_formats = c("analyzed_format")
  )
  permanova <- list()

  # Convert taxmap object to a phyloseq object.
  phyloseq_object <- metacoder::as_phyloseq(metacoder_object, otu_table = "otu_abundance", phy_tree = "phy_tree")

  # TODO: Convert this to metacoder strictly.
  permanova <- list()
  otu <- microbiome::abundances(phyloseq_object)
  meta <- microbiome::meta(phyloseq_object)

  # Use phyloseq to generate the distance if it is wunifrac/unifrac.
  if (distance_method == "wunifrac" | distance_method == "unifrac") {
    dist <- phyloseq::distance(phyloseq_object, method = distance_method)
    dist_formula <- as.formula(paste0("dist ~ ", group))
    permanova[["permanova"]] <- vegan::adonis(dist_formula, data = meta, permutations = 99)
  } else {
    dist_formula <- as.formula(paste0("t(otu) ~ ", group))
    permanova[["permanova"]] <- vegan::adonis(dist_formula, data = meta, permutations = 99, method = distance_method)
  }

  # Checking the homogeneity condition
  # Note the assumption of similar multivariate spread among the groups
  # ie. analogous to variance homogeneity
  # Here the groups have signif. different spreads and
  # permanova result may be potentially explained by that.
  dist <- vegan::vegdist(t(otu))
  permanova[["anova"]] <- anova(vegan::betadisper(dist, meta[[sym(group)]]))

  # Investigate the top factors
  # Show coefficients for the top taxa separating the groups
  permanova[["coefficients"]] <- coefficients(permanova$permanova)[paste0(group, "1"), ]
  if (is.null(permanova$coefficients)) {
    warning("Coefficients were not able to be generated using this distance method.")
  } else {
    permanova[["top_coefficients"]] <- permanova$coefficients[rev(order(abs(permanova$coefficients)))[1:50]]
  }
  return(permanova)
}

#' @title Top Coefficients Barplot
#' @description A barplot to visualize which taxa have had the greatest impact on the overall gut composition.
#' @param top_coefficients Coefficients for the top taxa separating the groups.
#' @param title The title of the plot, Default = NULL
#' @return Returns a barplot of the top coefficients
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(MicrobiomeR)
#'   data <- analyzed_silva
#'   p <- permanova(obj = data)
#'   top_coefficients_barplot(p$top_coefficients)
#' }
#' }
#' @export
#' @seealso View \code{\link{permanova}} to understand how to generate `top_coefficents`
#' @family Visualizations
#' @rdname top_coefficients_barplot
top_coefficients_barplot <- function(top_coefficients, title = NULL) {
  # Set graphical parameters
  par(mar = c(3, 16, 2, 2))

  if (is.null(title)) {
    title <- "Top taxa"
  }

  # Plot the top coefficients
  plot <- barplot(sort(top_coefficients),
    horiz = T, las = 1, main = title,
    col = ifelse(sort(top_coefficients) >= 0, "#3288bd", "#d53e4f")
  )
  return(plot)
}

#' @title Save Top Coefficients Barplot
#' @description Save a top coefficients barplot.
#' @param plot The plot object.
#' @param filename The name of the file. (an extension should not be included)
#' @param format The file format, Default = "tiff"
#' @param start_path The path of which your output should be saved.
#' @param ... Additional paramaters.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(MicrobiomeR)
#'   save_top_coefficients_barplot(plot = plot, filename = "stacked_bar_phylum")
#' }
#' }
#' @export
#' @rdname save_top_coefficients_barplot
#' @importFrom ggplot2 ggsave
#' @importFrom crayon green yellow
#' @importFrom glue glue
save_top_coefficients_barplot <- function(plot, filename, format = "tiff", start_path = "output", ...) {
  if (is.null(filename)) {
    filename <- plot
  }
  full_path <- output_dir(start_path = start_path, plot_type = "top_coef_barplot", ...)
  message(glue::glue(crayon::yellow("Saving Top Coefficents Barplot to the following directory: \n", "\r\t{full_path}")))
  message(glue::glue(crayon::green("Saving the top coefficients barplot.")))
  ggplot2::ggsave(paste0(filename, ".", format),
    plot = plot, device = format, path = full_path,
    width = 8, height = 5, units = "in", dpi = 500
  )
}
