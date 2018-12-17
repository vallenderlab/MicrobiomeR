## Initialize a list for storing heat tree plots
get_heat_tree_plots <- function(obj, rank_list = NULL, ...) {
  rank_index <- pkg.private$rank_index
  if (is.null(rank_list)) {
    rank_list <- unlist(pkg.private$ranks)
  }
  htrees <- list()
  # Create a metacoder object from a phyloseq/metacoder/RData file
  metacoder_object <- MicrobiomeR::object_handler(obj)
  metacoder_object <- MicrobiomeR::validate_MicrobiomeR_format(
    obj = metacoder_object,
    valid_formats = c("analyzed_format"),
    force_format = TRUE)
  # Create a list of heat_tree plots for saving
  for (rank in rank_list) {
    rank_level <- rank_index[[rank]]
    filtered_obj <- metacoder_object %>% taxa::filter_obs(data = c("statistical_data", "taxa_proportions", "taxa_abundance", "stats_tax_data"), n_supertaxa < rank_level, drop_taxa = TRUE)
    title <- sprintf("Bacterial Abundance (%s Level)", rank)
    default_heat_tree_parameters <- get_heat_tree_parameters(obj = filtered_obj, title = title, ...)
    # Filter by Taxonomy Rank and then create a heat tree.
    htrees[[rank]] <- do.call(what = "heat_tree", args = c(default_heat_tree_parameters))
    # Made plot title centered
    htrees[[rank]] <- htrees[[rank]] +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        text = ggplot2::element_text(size = 24, family = "Arial")) +
      ggplot2::ggtitle(title)
  }
  htrees[["metacoder_object"]] <- metacoder_object
  return(htrees)
}

get_heat_tree_parameters <- function(obj, title, ...) {
  input <- .data
  default_parameters <- list(
    .input = input,
    title = title,
    # NODE
    ## The node size is relevant to the Abundance level
    ## The node color is relevant to wheather the abundance is higher in control vs stressed animals
    ## The node labels are relevant to significant taxon names.
    node_size = input$n_obs(),
    node_color = log2_mean_ratio,
    node_label = ifelse(wilcox_p_value < 0.05, input$taxon_names(), NA),
    node_label_size = 1,
    ### The color red indicates higher abundance in Stressed animals
    ### The color blue indicates higher abundance in Control animals
    ### The color grey represents no difference in Control vs Stressed abundance
    ### Colorblind node_color_range = c("navy", "grey80", "greenyellow"),
    ### Colorblind node_color_range = c("navy", "grey80", "yellow3"),
    # #ffffbf
    node_color_range = c("#3288bd", "#f1f1f1", "#d53e4f"),
    node_color_trans = "linear",
    node_color_interval = c(-4, 4),
    node_size_axis_label = "Size: Number of OTUs",
    node_color_axis_label = "Color: Increased in Stressed (red) vs.\n Increased in Control (blue)\n",
    ### The labels are only for significant (pvalue < 0.05) abundance changes
    ### The labels are green to offset the blue/red colors.
    node_label_color = wilcox_p_value,
    node_label_color_range = c("darkgreen"),
    node_label_max = 1000,
    # EDGE (Branch)
    ## The edge size is relevant to the level of abundance.
    ## The edge color is relevant to significant levels of abundance.
    edge_color = wilcox_p_value,
    ### The color interval is set from 0 to 1 (p-value is a percentage).
    ### The color range is a vector consisting of a dark color followed
    ### by 19 of the same color.  This divides the interval by 20 so that
    ### every edge below 0.05 (5%; 1/20) is the dark color.
    ### So our p-value is the dark color.
    edge_color_range = c("honeydew2", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50", "grey50"),
    edge_color_trans = "linear",
    edge_color_interval = c(0, 1),
    edge_size_axis_label = "Size: Number of OTUs",
    edge_color_axis_label = "Color: Significant Changes \nAmong Treatments",
    # PLOT Options
    initial_layout = "reingold-tilford",
    layout = "davidson-harel",
    repel_labels = TRUE,
    repel_force = 3,
    overlap_avoidance = 3,
    make_edge_legend = FALSE
  )
  param_list <-list(params = default_parameters)
  param_list <- purrr::list_modify(param_list, ...)
  return(param_list)
}
