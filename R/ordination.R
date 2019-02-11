#' @title Ordination Plot
#' @description This function allows for ordination (which helps us to distinguish beta diversity relationships) to be plotted as well as for the corresponding data to be returned.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{create_metacoder}}.
#' @param method Choose an ordination method from 'PCoA', 'CCA', 'NMDS' or 'DPCoA', Default: 'PCoA'
#' @param distance Choose a distance method from 'bray', 'unifrac' or 'wunifrac', Default: 'wunifrac'
#' @param color Choose the group or factor of which colors will be mapped to, Default: 'TreatmentGroup'
#' @param title The title of the plot, Default: NULL
#' @param only_data Allows for only ordination data to be generated, Default: FALSE
#' @return By default, it returns an ordination plot.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # An example ordination plot
#'   library(MicrobiomeR)
#'   data <- analyzed_silva
#'   plot <- ordination_plot(obj = data)
#'   plot
#' }
#' }
#' @export
#' @family Visualizations
#' @rdname ordination_plot
#' @seealso View \code{\link{save_ordination_plots}} to save your ordination plot or multiple ordination plots.
#' @importFrom metacoder as_phyloseq
#' @importFrom phyloseq ordinate plot_ordination
#' @importFrom ggplot2 element_text geom_point theme element_blank guide_legend guides ggtitle unit scale_fill_manual labs scale_color_manual
ordination_plot <- function(obj, method = "PCoA", distance = "wunifrac", color = "TreatmentGroup", title = NULL, only_data = FALSE) {
  metacoder_object <- validate_MicrobiomeR_format(
    obj = create_metacoder(obj),
    valid_formats = c("analyzed_format")
  )

  # Convert metacoder object to a phyloseq object.
  phyloseq_object <- metacoder::as_phyloseq(metacoder_object, otu_table = "otu_abundance", phy_tree = "phy_tree")

  # Find out how many groups are in the "color" metadata
  meta <- phyloseq::sample_data(phyloseq_object)
  num_groups <- length(unique(meta[[color]]))
  group_names <- unique(meta[[color]])

  # Create palette and palette lists
  ord_palette <- c("#3288bd", "#d53e4f", "#62954C", "#C59144")
  names(ord_palette) <- group_names

  ord <- phyloseq::ordinate(phyloseq_object, method, distance)

  plot <- phyloseq::plot_ordination(physeq = phyloseq_object, ordination = ord, color = color) +
    ggplot2::geom_point(size = 2) + ggplot2::theme(
      text = ggplot2::element_text(size = 11, family = "Arial", face = "bold"),
      axis.text.y = ggplot2::element_text(margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
      axis.text.x = ggplot2::element_text(margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "cm")), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = element_blank(),
      axis.ticks.length = ggplot2::unit(-0.25, "cm"),
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.text = ggplot2::element_text(colour = "black"),
      panel.background = element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) + ggplot2::stat_ellipse(
      geom = "polygon", alpha = .2,
      ggplot2::aes(fill = !!sym(color)),
      show.legend = FALSE
    )

  plot <- plot + ggplot2::scale_fill_manual(values = ord_palette)

  plot <- plot + ggplot2::labs(color = color) + ggplot2::scale_color_manual(values = ord_palette)

  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title) + ggplot2::annotate("segment", x = Inf, xend = -Inf, y = Inf, yend = Inf, color = "black", lwd = 1) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(fill = NA)))
  }

  df <- phyloseq::plot_ordination(phyloseq_object, ord, justDF = TRUE, axes = 1:4)

  if (only_data == TRUE) {
    return(df)
  } else {
    return(plot)
  }
}

#' @title Ordination Plots
#' @description Generate plots for a list of ordination methods and distances.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{create_metacoder}}.
#' @param methods A list of ordination methods, Default: 'c("PCoA", "NMDS")'
#' @param distances A list of distance methods, Default: 'c("wunifrac", "unifrac", "bray")'
#' @return Returns a melted dataframe.
#' @family Visualizations
#' @rdname ordination_plot
ordination_plots <- function(obj, methods = c("PCoA", "NMDS"), distances = c("wunifrac", "unifrac", "bray"),
                             color = "TreatmentGroup", select_otu_table = "otu_proportions") {
  if (is.null(methods)) {
    methods <- c("PCoA", "NMDS")
  } else if (length(methods) < 2) {
    stop("Use the ordination_plot function for generating a plot for 1 method or 1 distance.")
  }
  if (is.null(distances)) {
    distances <- c("wunifrac", "unifrac", "bray")
  } else if (length(methods) < 2 && length(distances) < 2) {
    stop("Use the ordination_plot function for generating a plot for 1 method and 1 distance.")
  }
  ordination_plots <- list()
  for (m in methods) {
    for (d in distances) {
      ordination_plots[[paste0(m, "_", d)]] <- ordination_plot(obj, method = m, distance = d, color = color, title = NULL, only_data = FALSE)
    }
  }

  return(ordination_plots)
}


#' @title Save Ordination Plots
#' @description This function saves ordination plots storred in a listlike object to an output folder.
#' @param ord An ordination plot list generated by ordination_plots.
#' @param format The format of the output image.  Default: 'tiff'
#' @param start_path The starting path of the output directory.  Default: 'output'
#' @param ... An optional list of parameters to use in the output_dir function.
#' @return An output directory that contains ordination plots.
#' @pretty_print TRUE
#' @details This function creates an appropriate output directory, where it saves publication ready
#' plots.
#' @examples
#' \dontrun{
#' if (interactive()) {
#' }
#' }
#' @export
#' @family Visualizations
#' @rdname save_ordination_plots
#' @seealso
#'  \code{\link[MicrobiomeR]{output_dir}}
#'
#'  \code{\link[ggplot2]{ggsave}}
#' @importFrom ggplot2 ggsave
#' @importFrom crayon green
#' @importFrom glue glue
save_ordination_plots <- function(ord_plots, format = "tiff", start_path = "output", ...) {
  # Create the relative path to the ordination plots.  By default the path will be <pwd>/output/<experiment>/ordination/<format(Sys.time(), "%Y-%m-%d_%s")>
  # With the parameters set the full path will be <pwd>/output/<experiment>/ordination/<extra_path>.
  full_path <- output_dir(start_path = start_path, plot_type = "ordination", ...)
  message(glue::glue(crayon::yellow("Saving Ordination plots to the following directory: \n", "\r\t{full_path}")))
  # Iterate the plot list and save them in the proper directory
  for (method in names(ord_plots)) {
    if (method != "metacoder_object") {
      message(crayon::green("Saving the {method} ordination plot."))
      ggplot2::ggsave(filename = sprintf("%s_ordination.%s", method, format), plot = ord_plots[[method]], device = format, path = full_path, dpi = 500)
    }
  }
}
