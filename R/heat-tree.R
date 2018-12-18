#' @title Get Heat Tree Plots
#' @description A function for getting multiple heat_tree plots per rank.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param rank_list A vector of ranks used to generate heat_trees.  Default: NULL
#' @param ... Any of the \code{\link[metacoder]{heat_tree}} parameters can be used to change the way the heat_tree
#' output is displayed.
#' @return A list of heat_tree plots.
#' @pretty_print TRUE
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(MicrobiomeR)
#' h_trees <- get_heat_tree_plots(MicrobiomeR::analyzed_silva, rank_list = c("Phylum", "Class"))
#' h_trees$Class
#'  }
#' }
#' @export
#' @family Visualizations
#' @rdname get_heat_tree_plots
#' @seealso
#'  \code{\link[metacoder]{heat_tree}}
#'
#'  \code{\link[MicrobiomeR]{object_handler}},  \code{\link[MicrobiomeR]{validate_MicrobiomeR_format}},  \code{\link[MicrobiomeR]{get_heat_tree_parameters}}
#'
#'  \code{\link[taxa]{filter_obs}}
#'
#'  \code{\link[crayon]{crayon}}
#'
#'  \code{\link[ggplot2]{theme}},  \code{\link[ggplot2:element]{margin}},  \code{\link[ggplot2]{labs}}
#' @importFrom taxa filter_obs
#' @importFrom crayon bgWhite red
#' @importFrom metacoder heat_tree
#' @importFrom ggplot2 theme element_text ggtitle
get_heat_tree_plots <- function(obj, rank_list = NULL, ...) {
  rank_index <- pkg.private$rank_index
  if (is.null(rank_list)) {
    rank_list <- unlist(pkg.private$ranks)
  }
  htrees <- list()
  # Create a metacoder object from a phyloseq/metacoder/RData file
  metacoder_object <- object_handler(obj)
  metacoder_object <- validate_MicrobiomeR_format(
    obj = metacoder_object,
    valid_formats = c("analyzed_format"),
    force_format = TRUE)
  # Create a list of heat_tree plots for saving
  for (rank in rank_list) {
    rank_level <- rank_index[[rank]]
    filtered_obj <- metacoder_object %>% taxa::filter_obs(data = c("statistical_data", "taxa_proportions", "taxa_abundance", "stats_tax_data"), n_supertaxa < rank_level, drop_taxa = TRUE)
    title <- sprintf("Bacterial Abundance (%s Level)", rank)
    message(sprintf("Generating a Heat Tree for %s", crayon::bgWhite(crayon::red(title))))
    default_heat_tree_parameters <- get_heat_tree_parameters(obj = filtered_obj, title = title, ...)
    # Filter by Taxonomy Rank and then create a heat tree.
    htrees[[rank]] <- do.call(what = metacoder::heat_tree, args = default_heat_tree_parameters)
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


#' @title Get Heat Tree Parameters
#' @description This function get's the parameters used for the get_heat_tree_plots function.
#' @param obj A metacoder object.
#' @param title The title used in the heat_tree plot.
#' @param ... Any of the heat tree parameters list below can be used to change the way the heat_tree
#' output is displayed.  However, this function acts as a default list of parameters.  The memebers
#' of the default list will be overridden by the dot parameters.
#' @return A list used with do.call and the metacoder::heat_tree function.
#' @pretty_print TRUE
#' @export
#' @family Visualizations
#' @rdname get_heat_tree_parameters
#' @seealso
#'  \code{\link[metacoder]{heat_tree}}
#'
#'  \code{\link[taxa]{n_obs}},\code{\link[taxa]{taxon_names}}
#'  \code{\link[purrr]{list_modify}}
#' @importFrom taxa n_obs taxon_names
#' @importFrom purrr list_modify
get_heat_tree_parameters <- function(obj, title, ...) {
  input <- obj
  log2_mean_ratio <- input$data$statistical_data$log2_mean_ratio
  wilcox_p_value <- input$data$statistical_data$wilcox_p_value
  default_parameters <- list(
    .input = input,
    title = title,
    # NODE
    ## The node size is relevant to the Abundance level
    ## The node color is relevant to wheather the abundance is higher in control vs stressed animals
    ## The node labels are relevant to significant taxon names.
    node_size = taxa::n_obs(input),
    node_color = log2_mean_ratio,
    node_label = ifelse(wilcox_p_value < 0.05, taxa::taxon_names(input), NA),
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
  return(param_list$params)
}


#' @title Save Heat Tree Plots
#' @description This function saves heat tree plots storred in a list object to an output folder.
#' @param htrees A named list of heat trees.
#' @param format The format of the output image.  Default: 'tiff'
#' @param start_path The starting path of the output directory.  Default: 'output'
#' @param ... An optional list of parameters to use in the get_output_dir function.
#' @return An output directory that contains heat tree plots.
#' @pretty_print TRUE
#' @details This function creates an appropriate output directory, where it saves publication ready
#' plots.
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(MicrobiomeR)
#' h_trees <- get_heat_tree_plots(MicrobiomeR::analyzed_silva, rank_list = c("Phylum", "Class"))
#' # Save to \emph{./output/heat_trees} folder.
#' save_heat_tree_plots(h_trees)
#'  }
#' }
#' @export
#' @family Visualizations
#' @rdname save_heat_tree_plots
#' @seealso
#'  \code{\link[MicrobiomeR]{get_output_dir}}
#'  \code{\link[ggplot2]{ggsave}}
#' @importFrom ggplot2 ggsave
save_heat_tree_plots <- function (htrees, format="tiff", start_path = "output", ...) {
  # Create the relative path to the heat_tree plots.  By default the path will be <pwd>/output/<experiment>/heat_trees/<format(Sys.time(), "%Y-%m-%d_%s")>
  # With the parameters set the full path will be <pwd>/output/<experiment>/heat_trees/<extra_path>.
  full_path <- get_output_dir(start_path = start_path, plot_type = "heat_trees", ...)
  # Iterate the heat_tree plot list and save them in the proper directory
  for (rank in names(htrees)) {
    if (rank != "metacoder_object"){
      ggplot2::ggsave(filename = sprintf("%s.heat_tree.%s", rank, format), plot = htrees[[rank]], device = format, path = full_path, dpi = 500)
    }
  }
}
