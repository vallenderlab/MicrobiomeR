#' @title Correlation Plot
#' @description Create a correlation plot from a metacoder/taxmap object.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param primary_rank The primary rank used to label the points.
#' @param secondary_rank The secondary rank used to color the points.  Can be an integer specifying
#' the number of supertaxon ranks above the primary rank or the name of a supertaxon rank.  Default: TRUE
#' @param wp_value The Wilcoxian P-Value used to represent significant points.  Default: 0.05
#' @param pal_func A palette function that returns grDevices::colorRampPalette.
#' @return A 1:1 correlation plot built with ggplot2.
#' @pretty_print TRUE
#' @details Correlation plots help to better explain the heat tree findings.
#' @examples
#' \dontrun{
#' if(interactive()){
#' # This example uses data that are no longer available in the MicrobiomeR package,
#' # however, they can be easily generated with \code{\link{MicrobiomeR}{as_analyzed_format}}.
#' library(MicrobiomeR)
#' analyzed_silva <- as_MicrobiomeR_format(MicrobiomeR::raw_silva_2, "analyzed_format")
#' correlation_plot(analyzed_silva, primary_rank = "Class", secondary_rank = "Phylum")
#'  }
#' }
#' @export
#' @family Visualizations
#' @rdname correlation_plot
#' @seealso
#'  \code{\link[MicrobiomeR]{object_handler}},  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}},  \code{\link[MicrobiomeR]{correlation_data}},  \code{\link[MicrobiomeR]{plot_limits}},  \code{\link[MicrobiomeR]{get_color_palette}}
#'
#'  \code{\link[ggplot2]{ggplot}},  \code{\link[ggplot2]{aes}},  \code{\link[ggplot2]{geom_polygon}},  \code{\link[ggplot2]{geom_point}},  \code{\link[ggplot2]{labs}},  \code{\link[ggplot2]{scale_continuous}},  \code{\link[ggplot2]{scale_manual}},  \code{\link[ggplot2]{guide_legend}},  \code{\link[ggplot2]{geom_abline}}
#'
#'  \code{\link[ggrepel:geom_text_repel]{geom_label_repel}}
#'
#'  \code{\link[forcats]{fct_reorder}}
#' @importFrom ggplot2 ggplot aes geom_polygon geom_point labs scale_y_log10 scale_x_log10 scale_shape_manual scale_fill_manual scale_color_manual guide_legend geom_abline
#' @importFrom forcats fct_reorder
#' @importFrom ggrepel geom_label_repel
#' @importFrom crayon yellow
correlation_plot <- function(obj, primary_rank, secondary_rank = TRUE,
                             wp_value = 0.05, pal_func = NULL) {
  metacoder_object <- object_handler(obj)
  metacoder_object <- validate_MicrobiomeR_format(obj = metacoder_object,
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
      message(crayon::yellow("Your secondary_rank should be a higher rank than the label_taxa.  Their values are being switched."))
      secondary_rank <- primary_rank
      primary_rank <- secondary_rank
    }
  } else if (secondary_rank == FALSE) {
    secondary_rank <- primary_rank
  }

  # Get the corelation data
  corr_data <- correlation_data(obj = metacoder_object, primary_rank = primary_rank, secondary_rank = secondary_rank, wp_value = wp_value)
  corrs <- list()
  for (comp_title in names(corr_data)) {
    treat_data <- corr_data[[comp_title]]
    primary_data <- treat_data$primary_data
    significant_data <- treat_data$significant_data
    treatments <- treat_data$treatments

    # Get the limits of the plot based on the data
    plot_limits <- plot_limits(primary_data$mean_treat1, primary_data$mean_treat2)
    # Create a dataframe for background color
    background_limits <- data.frame(id = c("1", "1", "1", "2", "2", "2"), x = c(0, Inf, 0, 0, Inf, Inf), y = c(0, Inf, Inf, 0, 0, Inf))
    # Get a color palette
    secondary_taxa <- length(unique(primary_data[[(secondary_rank)]]))

    if (is.null(pal_func)) {
      pal_func <- combination_palette(
        magma = list(palette = viridis::magma, args = list(n=500), range=450:500, rev=TRUE),
        inferno = list(palette = viridis::inferno, args = list(n=500), range=100:400, rev=TRUE),
        cividis = list(palette = viridis::cividis, args = list(n=500), range=100:200, rev=TRUE),
        viridis = list(palette = viridis::viridis, args = list(n=500), range=100:480))
    }
    myPal <- get_color_palette(
      pal_func = pal_func,
      color_no = secondary_taxa,
      display = FALSE)

    # Start ggplot2 workflow
    corr <- ggplot2::ggplot(primary_data, ggplot2::aes(x = mean_treat1, y = mean_treat2)) +
      ggplot2::geom_polygon(background_limits, mapping = ggplot2::aes(x = x, y = y, fill = id), alpha = 0.07, show.legend = FALSE) +
      ggplot2::geom_point(data = significant_data, ggplot2::aes(shape = Abundance), size = 3.2, color = "black", stroke = 2, show.legend = FALSE) +
      ggplot2::geom_point(data = primary_data, ggplot2::aes(shape = Abundance, color = forcats::fct_reorder(primary_data[[secondary_rank]], color_wilcox_p_value, min)), size = 2.5, stroke = 1.5) +
      ggrepel::geom_label_repel(
        mapping = ggplot2::aes(label = rank_label), size = 3, segment.size = 0.15, point.padding = 0.5, box.padding = 0.6, alpha = 0.65, force = 45, max.iter = 10000, min.segment.length = 0.1, seed = 2289,
        nudge_x = ifelse(primary_data$mean_treat1 < primary_data$mean_treat2, -2, 2.5), nudge_y = ifelse(primary_data$mean_treat2 < primary_data$mean_treat1, -1.9, 2)) +
      ggplot2::labs(title = glue::glue("{primary_rank} ({comp_title})"), x = glue::glue("Mean Abundance Before {treatments[1]}"), y = glue::glue("Mean Abundance After {treatments[1]}")) +
      ggplot2::scale_y_log10(limits = plot_limits) + ggplot2::scale_x_log10(limits = plot_limits) +
      ggplot2::scale_shape_manual(name = glue::glue("Abundance After {treatments[1]}:"), values = c("Significant Increase" = 16, "Significant Decrease" = 15, "Insignificant Change" = 4)) +
      ggplot2::scale_fill_manual(values = c("red", "blue", myPal), guide = FALSE) +
      ggplot2::scale_color_manual(values = c(myPal), name = sprintf("%s:", c(secondary_rank)), guide = ggplot2::guide_legend(ncol = 2)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    message(crayon::green(sprintf("Generating Correlation Plot comparing %s for %s with color based on %s", crayon::bgWhite(comp_title), crayon::bgWhite(primary_rank), crayon::bgWhite(secondary_rank))))
    corrs[[comp_title]] <- corr
  }
  return(corrs)
}

#' @title Get Multiple Correlation Plots
#' @description This function allows the user to create a list of multiple correlation plots.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param primary_ranks A vector of primary ranks used to label the points.
#' @param secondary_ranks  The secondary rank used to color the points.  Can be an integer specifying
#' the number of supertaxon ranks above the primary rank or the name of a supertaxon rank.  Default: TRUE
#' @param ... An optional list of parameters to use in the correlation_plot function.
#' @return A list object containing correlation plots.  A pairwise comparison returns a nested list.
#' @pretty_print TRUE
#' @details This function makes it easier to generate multiple correlation plots at once.
#' @examples
#' \dontrun{
#' if(interactive()){
#' # This example uses data that are no longer available in the MicrobiomeR package,
#' # however, they can be easily generated with \code{\link{MicrobiomeR}{as_analyzed_format}}.
#' library(MicrobiomeR)
#' analyzed_silva <- as_MicrobiomeR_format(MicrobiomeR::raw_silva_2, "analyzed_format")
#' corr_plots <- correlation_plots(analyzed_silva, primary_ranks = c("Phylum", "Class", "Order"),
#'                       secondary_ranks = c("Phylum", "Class", "Order", "Family", "Genus"))
#' # Show a plot
#' corr_plots$Class$Phylum
#'  }
#' }
#' @export
#' @family Visualizations
#' @rdname correlation_plots
#' @importFrom crayon green bgWhite
correlation_plots <- function(obj, primary_ranks, secondary_ranks = TRUE, ...) {
  corr <- list()
  params <- list(...)
  rank_index <- pkg.private$rank_index
  ranks <- pkg.private$ranks
  # Allow secondary_ranks to be TRUE/FALSE isntead of vector
  if (is.vector(primary_ranks) && length(secondary_ranks) == 1) {
    # Get the primary and secondary rank of interest
    if (secondary_ranks == TRUE) {
      secondary_ranks <- 1
    }
    if (is.numeric(secondary_ranks) || is.character(secondary_ranks)) {
      secondary_ranks <- rep(secondary_ranks, length(primary_ranks))
    }
  }
  # Get correlation plots
  for (i in 1:length(primary_ranks)) {
    pr <- primary_ranks[i]
    corr[[pr]] <- list()
    for (j in 1:length(unique(secondary_ranks))) {
      sr <- secondary_ranks[j]
      if (is.numeric(sr)) {
        sr <- as.character(ranks[rank_index[[pr]] - sr])
      } else if (sr == FALSE) {
        sr <- pr
      }
      if (rank_index[[pr]] < rank_index[[sr]]) {
        next()
      }
      message(glue::glue(crayon::green("Generating a Correlation Plot comparing ", crayon::bgWhite({pr}), " with ", crayon::bgWhite({sr}), ".")))
      corr[[pr]][[sr]] <- do.call(correlation_plot, c(list(obj = obj, primary_rank = pr, secondary_rank = sr), params))
    }
  }
  return(corr)
}

#' @title Get Correlation Plot Data
#' @description Get the correlation plot data comparing all of the treatment groups.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param primary_rank A primary rank used to label the points.
#' @param secondary_rank  The secondary rank used to color the points, Default: TRUE
#' @param wp_value The wilcoxon p-value cutoff/threshold, Default: 0.05
#' @return Data for the correlation plots.
#' @pretty_print TRUE
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Visualizations
#' @rdname correlation_data
#' @seealso
#'  \code{\link[MicrobiomeR]{agglomerate_metacoder}},  \code{\link[MicrobiomeR]{vlookup}}
#'
#'  \code{\link[dplyr]{tidyeval}},  \code{\link[dplyr]{mutate}},  \code{\link[dplyr]{filter}},  \code{\link[dplyr]{arrange}}
#'
#'  \code{\link[taxa]{filter_obs}}
#' @importFrom dplyr enquo mutate filter arrange
#' @importFrom utils combn
#' @importFrom taxa filter_obs
#' @importFrom glue glue
correlation_data <- function(obj, primary_rank, secondary_rank = TRUE, wp_value = 0.05) {
  # Quotes
  quoted_str <- dplyr::enquo(secondary_rank)
  # Create the primary metacoder object
  primary_mo <- agglomerate_metacoder(obj = obj, rank = primary_rank,
                                      validated = TRUE)
  # Create the secondary metacoder object
  secondary_mo <- agglomerate_metacoder(obj = obj, rank = secondary_rank,
                                        validated = TRUE)
  # Get the primary and secondary statistical-taxonomy data frame
  primary_data <- primary_mo$data$stats_tax_data
  secondary_data <- secondary_mo$data$stats_tax_data

  combinations <- treatment_matrix(obj = obj)
  # Compare treatments
  corr_data <- list()
  for (i in seq_len(ncol(combinations))) {
    # Get treatment information
    treat_a <- combinations[1,i]
    treat_b <- combinations[2,i]

    p_mo <- primary_mo %>%
      taxa::filter_obs("stats_tax_data",
                       (treatment_1 == treat_a &
                          treatment_2 == treat_b) |
                         (treatment_1 == treat_b &
                            treatment_2 == treat_a))
    p_d <- p_mo$data$stats_tax_data
    s_mo <- secondary_mo %>%
      taxa::filter_obs("stats_tax_data",
                       (treatment_1 == treat_a &
                          treatment_2 == treat_b) |
                         (treatment_1 == treat_b &
                            treatment_2 == treat_a))
    s_d <- s_mo$data$stats_tax_data
    # Add a P-Value column for color, a Significance column (TRUE/FALSE) for subsetting, an Abundance column for the
    # legend, and a label column (TRUE/FALSE) for subsetting
    p_d <- p_d %>% dplyr::mutate(color_wilcox_p_value = vlookup(!! quoted_str, s_d, secondary_rank, "wilcox_p_value")) %>%
      dplyr::mutate(Significance = wilcox_p_value < wp_value) %>%
      dplyr::mutate(Abundance = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, "Significant Increase", "Significant Decrease"), "Insignificant Change")) %>%
      dplyr::mutate(label = ifelse(Significance == TRUE, ifelse(sign(log2_mean_ratio) == 1, TRUE, TRUE), FALSE))

    # Update points with 0 values so they look good on the plot
    p_d$mean_treat1[p_d$mean_treat1 == 0] <- 0.000000001
    p_d$mean_treat2[p_d$mean_treat2 == 0] <- 0.000000001

    # Create a label for the legend
    rank_label <- sprintf("%s_label", primary_rank)

    # Create a dataframe that contains only significant data.  This is used for point outlines
    sig_d <- dplyr::filter(p_d, label == TRUE) %>%
      dplyr::mutate(rank_label = .[[primary_rank]])

    # Add a column for labeling the points with taxonomy data based on significance
    p_d <- dplyr::mutate(p_d, rank_label = vlookup(p_d[[primary_rank]], sig_d, primary_rank, "rank_label"))
    p_d <- dplyr::arrange(p_d, label, wilcox_p_value)
    sig_d <- dplyr::arrange(sig_d, wilcox_p_value)
    corr_data[[glue::glue("{treat_a}-vs-{treat_b}")]] <- list(primary_data = p_d, significant_data = sig_d, treatments = combinations[,i])
  }
  return(corr_data)
}

#' @title Save Correlation Plots
#' @description This function saves correlation plots storred in a listlike object to an output folder.
#' @param corr A correlation plot list generated by correlation_plot or correlation_plots.
#' @param format The format of the output image.  Default: 'tiff'
#' @param start_path The starting path of the output directory.  Default: 'output'
#' @param ... An optional list of parameters to use in the output_dir function.
#' @return An output directory that contains correlation plots.
#' @pretty_print TRUE
#' @details This function creates an appropriate output directory, where it saves publication ready
#' plots.
#' @examples
#' \dontrun{
#' if(interactive()){
#' # This example uses data that are no longer available in the MicrobiomeR package,
#' # however, they can be easily generated with \code{\link{MicrobiomeR}{as_analyzed_silva}}.
#' library(MicrobiomeR)
#' analyzed_silva <- as_MicrobiomeR_format(MicrobiomeR::raw_silva_2, "analyzed_format")
#' corr_plot <- correlation_plot(analyzed_silva, primary_rank = "Class", secondary_rank = "Phylum")
#' # Save to \emph{./output/corr_plot} folder.
#' save_correlation_plots(corr_plot)
#'  }
#' }
#' @export
#' @family Visualizations
#' @rdname save_correlation_plots
#' @seealso
#'  \code{\link[MicrobiomeR]{output_dir}}
#'
#'  \code{\link[ggplot2]{ggsave}}
#' @importFrom ggplot2 ggsave
#' @importFrom crayon green
#' @importFrom glue glue
save_correlation_plots <- function(corr, format = "tiff", start_path = "output", ...) {
  # Create the relative path to the correlation plots.  By default the path will be <pwd>/output/<experiment>/corr_plot/<format(Sys.time(), "%Y-%m-%d_%s")>
  # With the parameters set the full path will be <pwd>/output/<experiment>/corr_plot/<extra_path>.
  full_path <- output_dir(start_path = start_path, plot_type = "corr_plot", ...)
  message(glue::glue(crayon::yellow("Saving Correlation Plots to the following directory: \n", "\r\t{full_path}")))
  # Iterate the heat_tree plot list and save them in the proper directory
  for (pr_name in names(corr)) {
    for (sr_name in names(corr[[pr_name]])) {
      for (comp in names(corr[[pr_name]][[sr_name]])) {
        message(glue::glue(crayon::green("Saving the {pr_name}-{sr_name} Correlation Plot comparing {comp}.")))
        ggplot2::ggsave(filename = sprintf("%s_%s_%s.corr_plot.tiff", pr_name, sr_name, comp), plot = corr[[pr_name]][[sr_name]][[comp]], device = format, path = full_path, dpi = 500, width = 500, height = 250, units = "mm")
      }
    }
  }
}


#' @title Get Plot Limits
#' @description Get the limits of the plot using data values.
#' @param x Data used to plot the x axis.
#' @param y Data used to plot the y axis.
#' @return A vector that supplies the overarching x and y limits.
#' @pretty_print TRUE
#' @family Visualizations
#' @rdname plot_limits
plot_limits <- function(x, y) {
  x_max <- max(x)
  x_min <- min(x)
  y_max <- max(y)
  y_min <- min(y)
  xy_limits <- c((min(c(x_min, y_min))), (max(c(x_max, y_max))))
  return(xy_limits)
}
