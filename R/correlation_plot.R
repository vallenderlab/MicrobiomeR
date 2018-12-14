correlation_plot <- function(obj, primary_rank, secondary_rank = TRUE,
                             wp_value = 0.05) {
  metacoder_object <- MicrobiomeR::object_handler(obj)
  metacoder_object <- MicrobiomeR::validate_MicrobiomeR_format(metacoder_object = metacoder_object,
                                                  valid_formats = c("analyzed_format"))
  ranks <- pkg.private$ranks
  rank_index <- pkg.private$rank_index
  # Get the parenk rank of the taxa of interest
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
  # Create the target data table
  primary_mo <- agglomerate_metacoder(metacoder_object = metacoder_object, rank = primary_rank,
                                      validated = TRUE)
  secondary_mo <- agglomerate_metacoder(metacoder_object = metacoder_object, rank = secondary_rank,
                                        validated = TRUE)
  target_data <- primary_mo$data$stats_tax_data
  super_taxa_data <- secondary_mo$data$stats_tax_data
  target_data <- target_data %>% mutate(color_wilcox_p_value = vlookup(!! quoted_str, super_taxa_data, secondary_rank, "wilcox_p_value")) %>%
    mutate(Significance = wilcox_p_value < wp_value) %>%
    mutate(Abundance = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, "Significant Increase", "Significant Decrease"), "Insignificant Change")) %>%
    mutate(label = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, TRUE, TRUE), FALSE))
  target_data$mean_treat1[target_data$mean_treat1 == 0] <- 0.000000001
  target_data$mean_treat2[target_data$mean_treat2 == 0] <- 0.000000001
  # Get perifieral data for the correlation plots
  limits <- get_plot_limits(target_data$mean_treat1, target_data$mean_treat2)
  positions <- data.frame(id = c("1", "1", "1", "2", "2", "2"), x = c(0, Inf, 0, 0, Inf, Inf), y = c(0, Inf, Inf, 0, 0, Inf))
  getPal <- colorRampPalette(brewer.pal(9, "Dark2"))
  myPal <- getPal(length(unique(target_data[[(secondary_rank)]])))
  # Get a dataframe for creating borders on the significant data
  rank_label <- sprintf("%s_label", primary_rank)
  ss_df <- dplyr::filter(target_data, label == TRUE) %>%
    dplyr::mutate(rank_label = .[[primary_rank]])
  target_data <- dplyr::mutate(target_data, rank_label = vlookup(target_data[[primary_rank]], ss_df, primary_rank, "rank_label"))
  target_data <- dplyr::arrange(target_data, label, wilcox_p_value)
  ss_df <- dplyr::arrange(ss_df, wilcox_p_value)

  # Start ggplot2 workflow
  corr <- ggplot(target_data, aes(x = mean_treat1, y = mean_treat2)) +
    geom_polygon(positions, mapping = aes(x = x, y = y, fill = id), alpha = 0.07, show.legend = FALSE) +
    geom_point(data = ss_df, aes(shape = Abundance), size = 3.2, color = "black", stroke = 2, show.legend = FALSE) +
    geom_point(data = target_data, aes(shape = Abundance, color = fct_reorder(target_data[[secondary_rank]], color_wilcox_p_value, min)), size = 2.5, stroke = 1.5) +
    geom_label_repel(
      mapping = aes(label = rank_label), size = 3, segment.size = 0.15, point.padding = 0.5, box.padding = 0.6, alpha = 0.65, force = 45, max.iter = 10000, min.segment.length = 0.1, seed = 2289,
      nudge_x = ifelse(target_data$mean_treat1 < target_data$mean_treat2, -2, 2.5), nudge_y = ifelse(target_data$mean_treat2 < target_data$mean_treat1, -1.9, 2)
    ) +
    labs(title = primary_rank, x = "Mean Abundance Before Stress", y = "Mean Abundance After Stress") +
    scale_y_log10(limits = limits) + scale_x_log10(limits = limits) +
    scale_shape_manual(name = "Abundance After Stress:", values = c("Significant Increase" = 16, "Significant Decrease" = 15, "Insignificant Change" = 4)) +
    scale_fill_manual(values = c("red", "blue", myPal), guide = FALSE) +
    scale_color_manual(values = c(myPal), name = sprintf("%s:", c(secondary_rank)), guide = guide_legend(ncol = 2)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  return(corr)

}
