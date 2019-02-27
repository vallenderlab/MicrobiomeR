#' @title Melt Taxmap
#' @param obj An object to be converted to a taxmap object with \code{\link[MicrobiomeR]{create_taxmap}}.
#' @importFrom dplyr right_join setdiff
#' @importFrom tidyr gather_
#' @family Formatting
#' @rdname melt_taxmap
melt_taxmap <- function(obj) {
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
#' @title Convert Proportions
#' @param melted_df A "melted" dataframe from the metacoder object's data.
#' @param tax_level The taxonomic level, Default: 'Phylum'
#' @importFrom dplyr filter group_by summarize mutate enquo quo_name
#' @importFrom stats na.omit
#' @family Data Manipulators
#' @rdname convert_proportions
convert_proportions <- function(melted_df, tax_level) {
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
#' @description Create a stacked barplot to show relative abundance of taxa. `convert_proportions` converts the dataframe abundance values to percent 100 and returns a transformed dataframe.
#' `melt_metacoder` melts the metacoder or phyloseq tables into a dataframe and returns a melted dataframe. `stacked_barplots` creates a stacked barplots for multiple taxonomic levels and returns a list of stacked barplots.
#' @param obj An object to be converted to a taxmap object with \code{\link[MicrobiomeR]{create_taxmap}}.
#' @param tax_level The taxonomic level, Default: 'Phylum'
#' @param fill The taxonomic level by which the bars are filled, Default: 'Phylum'
#' @param xlabel The label of the x axis, Default: 'Samples'
#' @param faceted A boolean to determine if the barplot should be faceted by TreatmentGroup
#' @param title The title or name of the plot.
#' @param palette_values A list of the colors to input to be mapped to the plot palette, Default: 'NULL'
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
#' @importFrom scales pretty_breaks
#' @importFrom shades scalefac saturation
#' @import vegan
#'
#' @inheritParams convert_proportions
#' @family Visualizations
#' @rdname stacked_barplot
#' @return Returns a stacked barplot.
#' @export
stacked_barplot <- function(obj, tax_level = "Phylum", fill = "Phylum", xlabel = "Samples", faceted = FALSE, title = NULL, palette_values = NULL) {
  metacoder_object <- validate_MicrobiomeR_format(
    obj = create_taxmap(obj),
    valid_formats = c("analyzed_format")
  )

  # Start by melting the data in the "standard" way using psmelt.
  # Also, transform the abundance data to relative abundance
  mdf <- convert_proportions(melt_taxmap(metacoder_object), tax_level)
  mdf <- dplyr::mutate(mdf, !!sym(tax_level) := factor(!!sym(tax_level), levels = unique(mdf[[tax_level]])))

  # Build the plot data structure
  p <- ggplot2::ggplot(mdf, aes(x = SampleID, y = Relative.Abundance, fill = !!sym(fill)), fill = fill)

  # Add the bar geometric object. Creates a basic graphic. Basis for the rest.
  p <- p + ggplot2::geom_bar(stat = "identity", position = "stack")

  # Add a palette if default not given.
  # Dynamically change palette colors based on number of taxa being input.
  if (is.null(palette_values)) {
    pal_func <- combination_palette(
      magma = list(palette = viridis::magma, args = list(n=500), range=450:500, rev=TRUE),
      inferno = list(palette = viridis::inferno, args = list(n=500), range=100:400, rev=TRUE),
      cividis = list(palette = viridis::cividis, args = list(n=500), range=100:200, rev=TRUE),
      viridis = list(palette = viridis::viridis, args = list(n=500), range=100:480))
    palette_values <- shades::saturation(get_color_palette(pal_func = pal_func, color_no = length(unique(mdf[[fill]])), display = FALSE),
                                         shades::scalefac(.6))
  }

  # Create the theme
  p <- p + ggplot2::ylab("Relative Abundance (% 100)") + ggplot2::scale_fill_manual(values = palette_values) + ggplot2::theme(
    text = ggplot2::element_text(size = 11, face = "bold"),
    axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "white"), strip.text = ggplot2::element_text(colour = "black"), panel.background = ggplot2::element_blank()
  ) +
    ggplot2::xlab(xlabel) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 4))

  # Add faceting, if true
  if (faceted) {
    stop("Faceted plots have not been fully integrated yet.")
    # p <- p + ggplot2::facet_grid(~TreatmentGroup, scales = "free_x", space = "free")
  }

  # Add a title, if given
  if (!is.null(title)) {
    p <- p + ggplot2::annotate("segment", x = Inf, xend = -Inf, y = Inf, yend = Inf, color = "black", lwd = 1) + ggplot2::ggtitle(title)
  }

  return(p)
}

#' @title Stacked Barplots
#' @param tax_levels The taxonomic levels, Default: 'c("Phylum", "Class", "Order")'
#' @param obj An object to be converted to a taxmap object with \code{\link[MicrobiomeR]{create_taxmap}}.
#' @family Visualizations
#' @rdname stacked_barplots
#' @export
stacked_barplots <- function(obj, tax_levels = c("Phylum", "Class", "Order")) {
  if (is.null(tax_levels)) {
    tax_levels <- c("Phylum", "Class", "Order")
  } else if (length(tax_levels) < 2) {
    stop("Use the stacked_plot function for generating a plot for 1 taxonomic level.")
  }
  stacked_barplots <- list()
  for (t in tax_levels) {
    stacked_barplots[[t]] <- stacked_barplot(obj = obj, tax_level = t, fill = t)
  }

  return(stacked_barplots)
}

#' @title Save Stacked Barplots
#' @description This function saves stacked barplot stored in a list object to an output folder.
#' @param sb_plots A named list of stacked barplots.
#' @param format The format of the output image.  Default: 'tiff'
#' @param start_path The starting path of the output directory.  Default: 'output'
#' @param ... An optional list of parameters to use in the output_dir function.
#' @return An output directory that contains stacked barplot.
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
#'   save_alpha_diversity_plot(salpha_div_plots)
#' }
#' }
#' @export
#' @family Visualizations
#' @rdname save_stacked_barplots
#' @seealso
#'  \code{\link[MicrobiomeR]{output_dir}}
#'
#'  \code{\link[ggplot2]{ggsave}}
#' @importFrom ggplot2 ggsave
#' @importFrom crayon yellow green
#' @importFrom glue glue
save_stacked_barplots <- function(sb_plots, format = "tiff", start_path = "output", ...) {
  # Create the relative path to the stacked barplots.  By default the path will be <pwd>/output/<experiment>/stacked_barplots/<format(Sys.time(), "%Y-%m-%d_%s")>
  # With the parameters set the full path will be <pwd>/output/<experiment>/stacked_barplots/<extra_path>.
  full_path <- output_dir(start_path = start_path, plot_type = "stacked_barplots", ...)
  message(glue::glue(crayon::yellow("Saving Stacked Barplots to the following directory: \n", "\r\t{full_path}")))
  # Iterate the plot list and save them in the proper directory
  for (rank in names(sb_plots)) {
    if (rank != "metacoder_object") {
      message(crayon::green("Saving the {rank} stacked barplot."))
      ggplot2::ggsave(
        filename = sprintf("%s_stacked_barplot.%s", tolower(rank), format), plot = sb_plots[[rank]], device = format, path = full_path,
        width = 8, height = 5, units = "in", dpi = 500
      )
    }
  }
}
