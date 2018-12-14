correlation_plot <- function(obj, primary_rank, secondary_rank = TRUE,
                             wp_value = 0.05) {
  metacoder_object <- MicrobiomeR::object_handler(obj)
  metacoder_object <- MicrobiomeR::validate_MicrobiomeR_format(metacoder_object = metacoder_object,
                                                  valid_formats = c("analyzed_format"))
  ranks <- pkg.private$ranks
  rank_index <- pkg.private$rank_index
  # Get the primary and secondary rank of interest
  if (secondary_rank == TRUE) {
    secondary_rank = 1
  }
  if (is.numeric(secondary_rank)) {
    secondary_rank <- as.character(ranks[rank_index[[primary_rank]] - secondary_rank])
  } else if (is.character(secondary_rank)) {
    if (rank_index[[secondary_rank]] < rank_index[[primary_rank]]) {
      secondary_rank <- secondary_rank
    } else {
      warning("Your secondary_rank should be a higher rank than the label_taxa.  Their values are being switched.")
      secondary_rank <- primary_rank
      primary_rank <- secondary_rank
    }
  } else if (secondary_rank == FALSE) {
    secondary_rank <- primary_rank
  }
  # Quotes
  quoted_str <- dplyr::enquo(secondary_rank)
  # Create the primary metacoder object
  primary_mo <- MicrobiomeR::agglomerate_metacoder(metacoder_object = metacoder_object, rank = primary_rank,
                                      validated = TRUE)
  # Create the secondary metacoder object
  secondary_mo <- MicrobiomeR::agglomerate_metacoder(metacoder_object = metacoder_object, rank = secondary_rank,
                                        validated = TRUE)
  # Get the primary and secondary statistical-taxonomy data frame
  primary_data <- primary_mo$data$stats_tax_data
  secondary_data <- secondary_mo$data$stats_tax_data

  # Add a P-Value column for color, a Significance column (TRUE/FALSE) for subsetting, an Abundance column for the
  # legend, and a label column (TRUE/FALSE) for subsetting
  primary_data <- primary_data %>% dplyr::mutate(color_wilcox_p_value = MicrobiomeR::vlookup(!! quoted_str, secondary_data, secondary_rank, "wilcox_p_value")) %>%
    dplyr::mutate(Significance = wilcox_p_value < wp_value) %>%
    dplyr::mutate(Abundance = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, "Significant Increase", "Significant Decrease"), "Insignificant Change")) %>%
    dplyr::mutate(label = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, TRUE, TRUE), FALSE))

  # Update points with 0 values so they look good on the plot
  primary_data$mean_treat1[primary_data$mean_treat1 == 0] <- 0.000000001
  primary_data$mean_treat2[primary_data$mean_treat2 == 0] <- 0.000000001

  # Create a label for the legend
  rank_label <- sprintf("%s_label", primary_rank)

  # Create a dataframe that contains only significant data.  This is used for point outlines
  significant_data <- dplyr::filter(primary_data, label == TRUE) %>%
    dplyr::mutate(rank_label = .[[primary_rank]])
  # Add a column for labeling the points with taxonomy data based on significance
  primary_data <- dplyr::mutate(primary_data, rank_label = MicrobiomeR::vlookup(primary_data[[primary_rank]], significant_data, primary_rank, "rank_label"))
  primary_data <- dplyr::arrange(primary_data, label, wilcox_p_value)
  significant_data <- dplyr::arrange(significant_data, wilcox_p_value)

  # Get the limits of the plot based on the data
  plot_limits <- MicrobiomeR::get_plot_limits(primary_data$mean_treat1, primary_data$mean_treat2)
  # Create a dataframe for background color
  background_limits <- data.frame(id = c("1", "1", "1", "2", "2", "2"), x = c(0, Inf, 0, 0, Inf, Inf), y = c(0, Inf, Inf, 0, 0, Inf))
  # Get a color palette
  getPal <- grDevices::colorRampPalette(brewer.pal(9, "Dark2"))
  myPal <- getPal(length(unique(primary_data[[(secondary_rank)]])))

  # Start ggplot2 workflow
  corr <- ggplot(primary_data, ggplot2::aes(x = mean_treat1, y = mean_treat2)) +
    ggplot2::geom_polygon(background_limits, mapping = ggplot2::aes(x = x, y = y, fill = id), alpha = 0.07, show.legend = FALSE) +
    ggplot2::geom_point(data = significant_data, ggplot2::aes(shape = Abundance), size = 3.2, color = "black", stroke = 2, show.legend = FALSE) +
    ggplot2::geom_point(data = primary_data, ggplot2::aes(shape = Abundance, color = forcats::fct_reorder(primary_data[[secondary_rank]], color_wilcox_p_value, min)), size = 2.5, stroke = 1.5) +
    ggrepel::geom_label_repel(
      mapping = ggplot2::aes(label = rank_label), size = 3, segment.size = 0.15, point.padding = 0.5, box.padding = 0.6, alpha = 0.65, force = 45, max.iter = 10000, min.segment.length = 0.1, seed = 2289,
      nudge_x = ifelse(primary_data$mean_treat1 < primary_data$mean_treat2, -2, 2.5), nudge_y = ifelse(primary_data$mean_treat2 < primary_data$mean_treat1, -1.9, 2)) +
    ggplot2::labs(title = primary_rank, x = "Mean Abundance Before Stress", y = "Mean Abundance After Stress") +
    ggplot2::scale_y_log10(limits = plot_limits) + ggplot2::scale_x_log10(limits = plot_limits) +
    ggplot2::scale_shape_manual(name = "Abundance After Stress:", values = c("Significant Increase" = 16, "Significant Decrease" = 15, "Insignificant Change" = 4)) +
    ggplot2::scale_fill_manual(values = c("red", "blue", myPal), guide = FALSE) +
    ggplot2::scale_color_manual(values = c(myPal), name = sprintf("%s:", c(secondary_rank)), guide = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  return(corr)

}
