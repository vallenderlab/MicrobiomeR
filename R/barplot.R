#' @title Melt Metacoder Object
#' @description Melt the metacoder or phyloseq tables into a dataframe.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @return Returns a melted dataframe.
#' @rdname melt_metacoder
#' @importFrom  dplyr right_join setdiff
#' @importFrom  tidyr gather_
melt_metacoder <- function(obj) {
  sd <- data.frame(obj$data$sample_data)
  TT <- data.frame(obj$data$otu_annotations, stringsAsFactors = FALSE)
  otu.table <- data.frame(obj$data$otu_proportions, check.names = FALSE, stringsAsFactors = FALSE)
  otu.table %>%
    dplyr::right_join(TT) %>%
    tidyr::gather("X.SampleID", "Abundance", dplyr::setdiff(colnames(otu.table), "otu_id")) %>%
    dplyr::right_join(sd) %>%
    rename(SampleID = `X.SampleID`) %>%
    rename(OTU = `otu_id`)
}

#' @title Transform Metacoder Dataframe
#' @description Transform the dataframe abundance values to percent 100.
#' @param melted_df A "melted" dataframe from the metacoder object's data.
#' @param tax_level The taxonomic level.
#' @return Returns a transformed dataframe.
#' @importFrom dplyr filter group_by summarize mutate enquo quo_name
#' @importFrom stats na.omit
transform_metacoder <- function(melted_df, tax_level) {
  t <- dplyr::enquo(tax_level)
  tax_level.abund <- paste0(dplyr::quo_name(t), ".Abundance")

  melted_df %>%
    dplyr::group_by(SampleID, !!sym(tax_level)) %>%
    dplyr::filter(Abundance > 0) %>%
    dplyr::summarize(!!tax_level.abund := sum(as.numeric(Abundance)), TreatmentGroup = dplyr::first(TreatmentGroup)) %>%
    stats::na.omit() %>%
    dplyr::mutate(Relative.Abundance = 100 * !!sym(tax_level.abund) / sum(!!sym(tax_level.abund)))
}

#' @title Stacked Barplot
#' @description Create a stacked barplot to show relative abundance of taxa.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param tax_level The taxonomic level, Default: 'Phylum'
#' @param fill The taxonomic level by which the bars are filled, Default: 'Phylum'
#' @param xlabel The label of the x axis, Default: 'Samples'
#' @param faceted A boolean to determine if the barplot should be faceted by TreatmentGroup
#' @param title The title or name of the plot.
#' @param palette_values A list of the colors to input to be mapped to the plot palette.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # An example stacked bar plot
#'   library(MicrobiomeR)
#'   data <- analyzed_silva
#'   palette <- get_color_palette(color_no = 12)
#'   plot <- stacked_barplot(obj = data, palette_values = palette)
#'   plot
#' }
#' }
#' @importFrom ggplot2 ggplot aes annotate geom_bar ylab element_blank element_rect xlab annotate
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @import scales
#' @import vegan
#'
#' @inheritParams transform_metacoder
#' @family Visualizations
#' @return Returns a stacked barplot.
#' @export
stacked_barplot <- function(obj, tax_level = "Phylum", fill = "Phylum", xlabel = "Samples", faceted = FALSE, title = NULL, palette_values) {
  metacoder_object <- validate_MicrobiomeR_format(
    obj = object_handler(obj),
    valid_formats = c("analyzed_format")
  )

  # Start by melting the data in the "standard" way using psmelt.
  # Also, transform the abundance data to relative abundance
  mdf <- transform_metacoder(melt_metacoder(metacoder_object), tax_level)
  mdf <- dplyr::mutate(mdf, !!sym(tax_level) := factor(!!sym(tax_level), levels = unique(mdf[[tax_level]])))

  # Build the plot data structure
  p <- ggplot2::ggplot(mdf, aes(x = SampleID, y = Relative.Abundance, fill = !!sym(fill)), fill = fill)

  # Add the bar geometric object. Creates a basic graphic. Basis for the rest.
  # Test weather additional
  p <- p + ggplot2::geom_bar(stat = "identity", position = "stack")

  # Add a title, if given
  if (!is.null(palette_values)) {
    # Dynamically change palette colors based on number of taxa being input.
    palette_values <- get_color_palette(display = FALSE)
  }

  # Create the theme
  p <- p + ggplot2::ylab("Relative Abundance (% 100)") + ggplot2::scale_fill_manual(values = palette_values) + ggplot2::theme(
    text = ggplot2::element_text(size = 11, face = "bold"),
    axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "white"), strip.text = ggplot2::element_text(colour = "black"), panel.background = ggplot2::element_blank()
  ) +
    ggplot2::xlab(xlabel) +
    ggplot2::annotate("segment", x = Inf, xend = -Inf, y = Inf, yend = Inf, color = "black", lwd = 1) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 4))

  # Add faceting, if true
  if (faceted) {
    stop("Faceted plots have not been fully integrated yet.")
    # p <- p + ggplot2::facet_grid(~TreatmentGroup, scales = "free_x", space = "free")
  }

  # Add a title, if given
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}

#' @title Save Barplot
#' @description Save a stacked barplot.
#' @param plot The plot object.
#' @param filename The name of the file. (an extension should not be included)
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(MicrobiomeR)
#'   save_barplot(plot = plot, filename = "stacked_bar_phylum")
#' }
#' }
#' @export
#' @rdname save_barplot
#' @importFrom ggplot2 ggsave
save_barplot <- function(plot, filename) {
  if (is.null(filename)) {
    filename <- plot
  }
  ggplot2::ggsave(paste0("output/", filename, ".tiff"),
    plot = plot, device = "tiff",
    width = 8, height = 5, units = "in", dpi = 500
  )
}
