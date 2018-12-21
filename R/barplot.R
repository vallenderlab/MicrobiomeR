#' @title Melt Metacoder Object
#' @description Melt the metacoder or phyloseq tables into a dataframe.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @importFrom dplyr right_join rename setdiff
#'
#' @return Returns a melted dataframe.
melt_metacoder_obj <- function(obj) {
  sd <- data.frame(obj$data$sample_data)
  TT <- data.frame(obj$data$otu_annotations, stringsAsFactors = FALSE)
  otu.table <- data.frame(obj$data$otu_proportions, check.names = FALSE, stringsAsFactors = FALSE)
  otu.table %>%
    right_join(TT) %>%
    gather_("X.SampleID", "Abundance", setdiff(colnames(otu.table), "otu_id")) %>%
    right_join(sd) %>%
    rename(SampleID = `X.SampleID`) %>%
    rename(OTU = `otu_id`)
}

#' @title Transform Metacoder Dataframe
#'
#' @description Transform the dataframe abundance values to percent 100.
#'
#' @param melted_df A "melted" dataframe from the metacoder object's data.
#' @param tax_level The taxonomic level.
#'
#' @importFrom dplyr filter group_by summarize mutate
#' @importFrom stats na.omit
#'
#' @return Returns a transformed dataframe.
transform_metacoder_df <- function(melted_df, tax_level) {
  # TODO: Add data wrangling step here or object validation.
  t <- enquo(tax_level)
  tax_level.abund <- paste0(quo_name(t), ".Abundance")

  melted_df %>%
    dplyr::group_by(SampleID, !!sym(tax_level)) %>%
    dplyr::filter(Abundance > 0) %>%
    dplyr::summarize(!!tax_level.abund := sum(as.numeric(Abundance)), TreatmentGroup = first(TreatmentGroup)) %>%
    stats::na.omit() %>%
    dplyr::mutate(Relative.Abundance = 100 * !!sym(tax_level.abund) / sum(!!sym(tax_level.abund)))
}

#' @title Stacked Barplot
#'
#' @description FUNCTION_DESCRIPTION
#'
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param tax_level DESCRIPTION.
#' @param fill DESCRIPTION.
#' @param xlabel DESCRIPTION.
#' @param title The title or name of the plot.
#' @param palette_values DESCRIPTION.
#'
#' @importFrom ggplot2 ggplot aes annotate geom_bar
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#'
#' @inheritParams transform_metacoder_obj
#' @family Visualizations
#' @return Returns a stacked barplot.
#' @export
stacked_barplot <- function(obj, tax_level = "Phylum", fill = "Phylum", xlabel = "Samples", title = NULL, palette_values) {

  #
  metacoder_object <- object_handler(obj)
  metacoder_object <- validate_MicrobiomeR_format(
    obj = metacoder_object,
    valid_formats = c("analyzed_format")
  )

  # Start by melting the data in the "standard" way using psmelt.
  # Also, transform the abundance data to relative abundance
  mdf <- transform_metacoder_df(melt_metacoder_obj(metacoder_object), tax_level)
  mdf <- dplyr::mutate(mdf, !!sym(tax_level) := factor(!!sym(tax_level), levels = unique(mdf[[tax_level]])))

  # Build the plot data structure
  p <- ggplot2::ggplot(mdf, aes(x = SampleID, y = Relative.Abundance, fill = !!sym(fill)), fill = fill)

  # Add the bar geometric object. Creates a basic graphic. Basis for the rest.
  # Test weather additional
  p <- p + ggplot2::geom_bar(stat = "identity", position = "stack")

  # Add a title, if given
  if (!is.null(palette_values)) {
    # Dynamically change palette colors based on number of taxa being input.
    palette_values <- MicrobiomeR::get_color_palette()
  }

  # Create the theme
  p <- p + ylab("Relative Abundance (% 100)") + scale_fill_manual(values = palette_values) + theme(
    text = element_text(size = 11, family = "Arial", face = "bold"),
    axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"), strip.text = element_text(colour = "black"), panel.background = element_blank()
  ) +
    ggplot2::xlab(xlabel) +
    ggplot2::annotate("segment", x = Inf, xend = -Inf, y = Inf, yend = Inf, color = "black", lwd = 1) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))

  # Add faceting
  p <- p + facet_grid(~TreatmentGroup, scales = "free_x", space = "free")

  # Add a title, if given
  if (!is.null(title)) {
    p <- p + ggtitle(title)
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
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname save_barplot
#' @seealso
#'
#' @import ggplot2
save_barplot <- function(plot, filename) {
  if (is.null(filename)) {
    filename <- plot
  }
  ggplot2::ggsave(paste0("output/", filename, ".tiff"),
    plot = plot, device = "tiff",
    width = 8, height = 5, units = "in", dpi = 500
  )
}
