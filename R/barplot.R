#' @title Melt Metacoder Object
#'
#' @description Melt the metacoder or phyloseq tables into a dataframe.
#'
#' @param obj A metacoder or phyloseq object.
#'
#' @importFrom dplyr right_join rename gather_ setdiff
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
#' @param melted_df DESCRIPTION.
#' @param tax_level DESCRIPTION.
#'
#' @importFrom dplyr group_by summarize filter mutate
#'
#' @return Returns a transformed dataframe.
transform_metacoder_df <- function(melted_df, tax_level) {
  # TODO: Add data wrangling step here or object validation.
  t <- enquo(tax_level)
  tax_level.abund <- paste0(quo_name(t), ".Abundance")

  melted_df %>%
    group_by(SampleID, !!sym(tax_level)) %>%
    filter(Abundance > 0) %>%
    summarize(!!tax_level.abund := sum(as.numeric(Abundance)), TreatmentGroup = first(TreatmentGroup)) %>%
    na.omit() %>%
    mutate(Relative.Abundance = 100 * !!sym(tax_level.abund) / sum(!!sym(tax_level.abund)))
}



#' @title Stacked Barplot
#'
#' @description FUNCTION_DESCRIPTION
#'
#' @param obj A metacoder or phyloseq object.
#' @param tax_level DESCRIPTION.
#' @param fill DESCRIPTION.
#' @param xlabel DESCRIPTION.
#' @param title DESCRIPTION.
#' @param palette_values DESCRIPTION.
#'
#' @importFrom ggplot2 ggplot aes annotate geom_bar
#' @importFrom magrittr %>%
#'
#' @inheritParams transform_metacoder_obj
#'
#' @return Returns a stacked barplot.
#' @export
stacked_barplot <- function(obj, tax_level = "Phylum", fill = "Phylum", xlabel = "Samples", title = NULL, palette_values) {

  # TODO: Add data wrangling step here or object validation.
  # TODO: Add default palette values.

  # Start by melting the data in the "standard" way using psmelt.
  # Also, transform the abundance data to relative abundance
  mdf <- transform_metacoder_df(melt_metacoder_obj(obj), tax_level)
  mdf <- mutate(mdf, !!sym(tax_level) := factor(!!sym(tax_level), levels = unique(mdf[[tax_level]])))

  # Build the plot data structure
  p <- ggplot2::ggplot(mdf, aes(x = SampleID, y = Relative.Abundance, fill = !!sym(fill)), fill = fill)

  # Add the bar geometric object. Creates a basic graphic. Basis for the rest.
  # Test weather additional
  p <- p + geom_bar(stat = "identity", position = "stack")

  # Create the theme
  p <- p + ylab("Relative Abundance (% 100)") + scale_fill_manual(values = palette_values) + theme(
    text = element_text(size = 11, family = "Arial", face = "bold"),
    axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"), strip.text = element_text(colour = "black"), panel.background = element_blank()
  ) +
    xlab(xlabel) +
    annotate("segment", x = Inf, xend = -Inf, y = Inf, yend = Inf, color = "black", lwd = 1) +
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
#'
#' @description Save a stacked barplot.
#'
#' @param plot The stacked barplot object.
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
