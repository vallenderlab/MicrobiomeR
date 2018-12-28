#' @title Ordination Plot
#' @description This function allows for ordination (which helps us to distinguish beta diversity relationships) to be plotted as well as for the corresponding data to be returned.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param method Choose an ordination method from 'PCoA', 'CCA', 'NMDS' or 'DPCoA', Default: 'PCoA'
#' @param distance Choose a distance method from 'bray', 'unifrac' or 'wunifrac', Default: 'wunifrac'
#' @param color Choose the group or factor of which colors will be mapped to, Default: 'TreatmentGroup'
#' @param title The title of the plot, Default: NULL
#' @param save Save your plot, Default: FALSE
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
#' @rdname ordination_plot
#' @seealso
#'
#' @importFrom metacoder as_phyloseq
#' @importFrom phyloseq ordinate plot_ordination
#' @importFrom ggplot2 element_text geom_point theme element_blank guide_legend guides ggtitle unit scale_fill_manual labs scale_color_manual
ordination_plot <- function(obj, method = "PCoA", distance = "wunifrac", color = "TreatmentGroup", title = NULL, save = FALSE, only_data = FALSE) {
  metacoder_object <- validate_MicrobiomeR_format(
    obj = object_handler(obj),
    valid_formats = c("analyzed_format")
  )

  # Convert metacoder object to a phyloseq object.
  phyloseq_object <- metacoder::as_phyloseq(metacoder_object, otu_table = "otu_abundance", phy_tree = "phy_tree")

  ord <- phyloseq::ordinate(phyloseq_object, method, distance)

  plot <- phyloseq::plot_ordination(physeq = phyloseq_object, ordination = ord, color = color) +
    ggplot2::geom_point(size = 2) + ggplot2::theme(
      text = ggplot2::element_text(size = 11, family = "Arial", face = "bold"),
      axis.text.y = ggplot2::element_text(margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
      axis.text.x = ggplot2::element_text(margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "cm")), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = element_blank(),
      axis.ticks.length = ggplot2::unit(-0.25, "cm"),
      strip.background = ggplot2::element_rect(fill = "white"), strip.text = ggplot2::element_text(colour = "black"), panel.background = element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) + ggplot2::stat_ellipse(geom = "polygon", alpha = .2, ggplot2::aes(fill = !!sym(color)), show.legend = FALSE)
      + ggplot2::scale_fill_manual(values = c("#3288bd", "#d53e4f")) # TODO: Use a palette here instead.

  # TODO: Add the ability to change variables here or use multiple variables. Could be better to remove this altogether.
  plot <- plot + ggplot2::labs(color = color) + ggplot2::scale_color_manual(values = c("Control" = "#3288bd", "Stressed" = "#d53e4f"))

  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title) + ggplot2::annotate("segment", x = Inf, xend = -Inf, y = Inf, yend = Inf, color = "black", lwd = 1) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(fill = NA)))
  }

  df <- phyloseq::plot_ordination(phyloseq_object, ord, justDF = TRUE, axes = 1:4)

  if (save == TRUE) {
    ggplot2::ggsave(filename = paste0("ordination_", tolower(method), "_", tolower(distance), ".tiff"))
    return(plot)
  } else if (only_data == TRUE) {
    return(df)
  } else if (save == FALSE && only_data == TRUE) {
    warning("You may only return the data or the plot. The plot will be returned by default.")
    return(plot)
  } else if (save == TRUE && only_data == TRUE) {
    warning("Since you are saving the plot, only the data will be returned by default.")
    return(df)
  } else {
    return(plot)
  }
}
