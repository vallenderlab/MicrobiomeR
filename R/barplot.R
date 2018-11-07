#' @title Facet Bar Plot
#'
#' @description A faceted stacked bar plot ideal for plotting 2 groups.
#'
#' @param phyloseq_object A phyloseq object.
#' @param fill
#' @param title
#' @param xlabel
#'
#' @return Returns a stacked barplot.
#' @export
facet_barplot <- function(phyloseq_object, fill, title, xlabel = "Samples") {
  phyloseq::plot_bar(phyloseq_object, fill = fill) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::ylab("Relative Abundance (% 100)") +
    ggplot2::facet_grid(cols = ggplot2::vars(TreatmentGroup), scales = "free_x", space = "free") +
    ggplot2::theme(
      text = ggplot2::element_text(size = 11, family = "Arial", face = "bold"),
      axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white"), strip.text = ggplot2::element_text(colour = "black"), panel.background = ggplot2::element_blank()
    ) +
    ggplot2::xlab(xlabel) + ggplot2::ggtitle(title) +
    ggplot2::annotate("segment", x = Inf, xend = -Inf, y = Inf, yend = Inf, color = "black", lwd = 1) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
}

#' @title Save Barplot
#'
#' @description Save a stacked barplot.
#'
#' @param plot
#' @param filename The filename of the plot. (It will be saved as .tiff)
#'
#' @export
save_barplot <- function(plot, filename) {
  if (is.null(filename)) {
    filename <- plot
  }
  ggplot2::ggsave(paste0("output/", filename, ".tiff"),
         plot = plot, device = "tiff",
         width = 8, height = 5, units = "in", dpi = 500
  )
}
