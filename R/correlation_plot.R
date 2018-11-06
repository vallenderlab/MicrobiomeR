
# #TODO: Rename correlation_plot().
# TODO: Create function for generating plotly data.
# TODO: Create function for generating shiny data table.
# TODO: Call functions with :: for package.
#' @title Get a Correlation Plot
#' @description FUNCTION_DESCRIPTION
#' @param phyloseq_object PARAM_DESCRIPTION
#' @param label_rank PARAM_DESCRIPTION
#' @param super_taxa PARAM_DESCRIPTION, Default: 1
#' @param wp_value PARAM_DESCRIPTION, Default: 0.05
#' @param plotly PARAM_DESCRIPTION, Default: FALSE
#' @param excel_export PARAM_DESCRIPTION, Default: FALSE
#' @param plotly_export PARAM_DESCRIPTION, Default: FALSE
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @pretty_print TRUE
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family FAMILY_TITLE
#' @rdname correlation_plot
#' @seealso
#'  \code{\link[phyloseq]{tax_glom}}
#'  \code{\link[dplyr]{mutate}},\code{\link[dplyr]{filter}},\code{\link[dplyr]{arrange}}
#'  \code{\link[stringr]{str_detect}}
#'  \code{\link[plotly]{ggplotly}}
#'  \code{\link[DT]{dataTableOutput}}
#' @importFrom phyloseq tax_glom
#' @importFrom dplyr mutate filter arrange
#' @importFrom stringr str_detect
#' @importFrom plotly ggplotly
#' @importFrom DT dataTableOutput renderDataTable
correlation_plot <- function(phyloseq_object, label_rank, super_taxa = 1, wp_value = 0.05,
                                 plotly = FALSE, excel_export = FALSE, plotly_export = FALSE, ...) {
  # TODO:  Create rank_index/rank variables dynamically or from data
  if (is.numeric(super_taxa)) {
    super_rank <- as.character(ranks[rank_index[[label_rank]] - super_taxa])
  } else if (is.character(super_taxa)) {
    if (rank_index[[super_taxa]] < rank_index[[label_rank]]) {
      super_rank <- super_taxa
    } else {
      warning("Your super_taxa should be a higher rank than the label_taxa.  Their values are being switched.")
      super_rank <- label_rank
      label_rank <- super_taxa
    }
  }
  # PHYLOSEQ data manipulation
  label_phy <- phyloseq::tax_glom(phyloseq_object, taxrank = label_rank)
  color_phy <- phyloseq::tax_glom(phyloseq_object, taxrank = super_rank)

  label_data <- phyloseq_to_stats_dataframe(label_phy, c("Stressed", "Control"))
  color_data <- phyloseq_to_stats_dataframe(color_phy, c("Stressed", "Control"))
  label_data$data <- label_data$data %>% mutate(color_wilcox_p_value = vlookup(Phylum, color_data$data, "Phylum", "wilcox_p_value"))
  # Create columns for easier plotting
  df_tax_otu <- label_data$data
  df_tax_otu <- mutate(df_tax_otu, Significance = wilcox_p_value < wp_value)
  df_tax_otu <- dplyr::mutate(df_tax_otu, Abundance = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, "Significant Increase", "Significant Decrease"), "Insignificant Change"))
  df_tax_otu <- dplyr::mutate(df_tax_otu, label = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, TRUE, TRUE), FALSE))
  if (excel_export == TRUE) {
    # TODO: Update the phylsoeq_to_excel function to take a dataframe and list of column names.
    # phyloseq_to_excel()
  }
  limits <- get_plot_limits(df_tax_otu$mean_stressed, df_tax_otu$mean_control)
  rank_label <- sprintf("%s_label", label_rank)
  # Further data manipulation
  ss_df <- dplyr::filter(df_tax_otu, label == TRUE)
  ss_df <- dplyr::filter(ss_df, !stringr::str_detect(ss_df[[label_rank]], "uncultured"))
  ss_df <- dplyr::mutate(ss_df, rank_label = ss_df[[label_rank]])
  ss_df <- dplyr::arrange(ss_df, wilcox_p_value)
  df_tax_otu <- dplyr::mutate(df_tax_otu, rank_label = vlookup(df_tax_otu[[label_rank]], ss_df, label_rank, "rank_label"))
  df_tax_otu <- dplyr::arrange(df_tax_otu, label, wilcox_p_value)
  # Color Palette creation and geom parameter creation
  getPal <- colorRampPalette(brewer.pal(9, "Dark2"))
  myPal <- getPal(length(unique(df_tax_otu[[(super_rank)]])))
  positions <- data.frame(id = c("1", "1", "1", "2", "2", "2"), x = c(0, Inf, 0, 0, Inf, Inf), y = c(0, Inf, Inf, 0, 0, Inf))

  # Plot
  # TODO: Remove the Control/Stress items from correlation plot.  Make dynamic.
  corr <- ggplot(df_tax_otu, aes(x = mean_control, y = mean_stressed))
  # Add background color to further distinguish treatment comparison
  corr <- corr + geom_polygon(positions, mapping = aes(x = x, y = y, fill = id), alpha = 0.07, show.legend = FALSE)
  # Add first point for highlighting significant taxonomies
  # Add second point for giving the point data color based on the super_rank
  corr <- corr + geom_point(data = ss_df, aes(shape = Abundance), size = 3.2, color = "black", stroke = 2, show.legend = FALSE) +
    geom_point(data = df_tax_otu, aes(shape = Abundance, color = fct_reorder(df_tax_otu[[super_rank]], color_wilcox_p_value, min)), size = 2.5, stroke = 1.5)
  # Add labels to the significant points and try to space them out properly
  # TODO: Use this link for managing labels https://mikewk.com/post/2018-09-20-labelling-dataviz/
  corr <- corr + geom_label_repel(mapping = aes(label = rank_label), size = 3, segment.size = 0.15, point.padding = 0.5, box.padding = 0.6, alpha = 0.65, force = 45, max.iter = 10000, min.segment.length = 0.1, seed = 2289,
                     nudge_x = ifelse(df_tax_otu$mean_control < df_tax_otu$mean_stressed, -2, 2.5), nudge_y = ifelse(df_tax_otu$mean_stressed < df_tax_otu$mean_control, -1.9, 2))
  # Add titles and scale x/y axis
  corr <- corr + labs(title = label_rank, x = "Mean Abundance Before Stress", y = "Mean Abundance After Stress") +
    scale_y_log10(limits = limits) + scale_x_log10(limits = limits)
  # Add shapes based on significance
  corr <- corr + scale_shape_manual(name = "Abundance After Stress:", values = c("Significant Increase" = 16, "Significant Decrease" = 15, "Insignificant Change" = 4))
  # Add color to the points?  and polygons
  # TODO: Test if the points can be removed.
  corr <- corr + scale_fill_manual(values = c("red", "blue", myPal), guide = FALSE)
  # Add color to the points based on super_rank
  corr <- corr + scale_color_manual(values = c(myPal), name = sprintf("%s:", c(super_rank)), guide=guide_legend(ncol = 2))
  # Add a 1:1 ratio line
  corr <- corr + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  return(corr)
}

point_data <- function(phyloseq_obj, label_rank, color_rank) {
  # PHYLOSEQ data manipulation
  label_phy <- phyloseq::tax_glom(phyloseq_object, taxrank = label_rank)
  color_phy <- phyloseq::tax_glom(phyloseq_object, taxrank = color_rank)

  label_data <- phyloseq_to_stats_dataframe(label_phy, c("Stressed", "Control"))
  color_data <- phyloseq_to_stats_dataframe(color_phy, c("Stressed", "Control"))
  label_data$data <- label_data$data %>% mutate(color_wilcox_p_value = vlookup(Phylum, color_data$data, "Phylum", "wilcox_p_value"))
  # Create columns for easier plotting
  df_tax_otu <- label_data$data
  df_tax_otu <- mutate(df_tax_otu, Significance = wilcox_p_value < wp_value)
  df_tax_otu <- dplyr::mutate(df_tax_otu, Abundance = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, "Significant Increase", "Significant Decrease"), "Insignificant Change"))
  df_tax_otu <- dplyr::mutate(df_tax_otu, label = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, TRUE, TRUE), FALSE))
  if (excel_export == TRUE) {
    # TODO: Update the phylsoeq_to_excel function to take a dataframe and list of column names.
    # phyloseq_to_excel()
  }
  limits <- get_plot_limits(df_tax_otu$mean_stressed, df_tax_otu$mean_control)
  rank_label <- sprintf("%s_label", label_rank)
  # Further data manipulation
  ss_df <- dplyr::filter(df_tax_otu, label == TRUE)
  ss_df <- dplyr::filter(ss_df, !stringr::str_detect(ss_df[[label_rank]], "uncultured"))
  ss_df <- dplyr::mutate(ss_df, rank_label = ss_df[[label_rank]])
  ss_df <- dplyr::arrange(ss_df, wilcox_p_value)
  df_tax_otu <- dplyr::mutate(df_tax_otu, rank_label = vlookup(df_tax_otu[[label_rank]], ss_df, label_rank, "rank_label"))
  df_tax_otu <- dplyr::arrange(df_tax_otu, label, wilcox_p_value)
  # Color Palette creation and geom parameter creation
  getPal <- colorRampPalette(brewer.pal(9, "Dark2"))
  myPal <- getPal(length(unique(df_tax_otu[[(super_rank)]])))
  positions <- data.frame(id = c("1", "1", "1", "2", "2", "2"), x = c(0, Inf, 0, 0, Inf, Inf), y = c(0, Inf, Inf, 0, 0, Inf))

}
