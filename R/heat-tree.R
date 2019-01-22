#' @title Get Heat Tree Plots
#' @description A function for getting multiple heat_tree plots per rank.
#' @param obj An object to be converted to a metacoder object with \code{\link[MicrobiomeR]{object_handler}}.
#' @param rank_list A vector of ranks used to generate heat_trees.  Default: NULL
#' @param ... Any of the \code{\link[metacoder]{heat_tree}} parameters can be used to change the way the heat_tree
#' output is displayed.  Please see the \code{\link[MicrobiomeR]{get_heat_tree_parameters}} documentation
#' for further explanation.
#' @return A list of heat_tree plots.
#' @pretty_print TRUE
#' @examples
#' \dontrun{
#' if(interactive()){
#' # This example uses data that are no longer available in the MicrobiomeR package,
#' # however, they can be easily generated with \code{\link{MicrobiomeR}{as_analyzed_format}}.
#' library(MicrobiomeR)
#' analyzed_silva <- as_MicrobiomeR_format(MicrobiomeR::raw_silva_2, "analyzed_format")
#' h_trees <- get_heat_tree_plots(analyzed_silva, rank_list = c("Phylum", "Class"))
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
#' @importFrom crayon green bgWhite
get_heat_tree_plots <- function(obj, rank_list = NULL, ...) {
  suppressWarnings({
    rank_index <- pkg.private$rank_index
    if (is.null(rank_list)) {
      rank_list <- unlist(pkg.private$ranks)
    }
    returned_data <- list()
    htrees <- list()
    flt_taxmaps <- list()
    # Create a metacoder object from a phyloseq/metacoder/RData file
    obj <- object_handler(obj)
    obj <- validate_MicrobiomeR_format(
      obj = obj,
      valid_formats = c("analyzed_format"),
      force_format = TRUE)
    returned_data[["metacoder_object"]] <- obj
    # Create a list of heat_tree plots for saving
    for (rank in rank_list) {
      rank_level <- rank_index[[rank]]
      filtered_obj <- obj %>% taxa::filter_taxa(n_supertaxa < rank_level,
                                                supertaxa = TRUE,
                                                reassign_obs = FALSE)
      flt_taxmaps[[rank]] <- filtered_obj
      title <- sprintf("Bacterial Abundance (%s Level)", rank)
      message(crayon::green(sprintf("Generating a Heat Tree for %s", crayon::bgWhite(title))))
      treatment_no <- length(unique(filtered_obj$data$sample_data$TreatmentGroup))
      default_heat_tree_parameters <- get_heat_tree_parameters(obj = filtered_obj, title = title, treatment_no = treatment_no, ...)
      # Filter by Taxonomy Rank and then create a heat tree.
      if (treatment_no == 2) {
        htrees[[rank]] <- do.call(what = metacoder::heat_tree, args = default_heat_tree_parameters)
      } else if (treatment_no > 2) {
        #return(default_heat_tree_parameters)
        htm <- function(...) {
          do.call(what = metacoder::heat_tree_matrix, args = c(list(obj = filtered_obj, data="statistical_data"), ...))
        }
        htrees[[rank]] <- htm(default_heat_tree_parameters)
      }
      # Made plot title centered
      htrees[[rank]] <- htrees[[rank]] +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          text = ggplot2::element_text(size = 24, family = "Arial"))
    }
    returned_data[["heat_trees"]] <- htrees
    returned_data[["taxmaps"]] <- flt_taxmaps
  })
  return(returned_data)
}


#' @title Get Heat Tree Parameters
#' @description This function get's the parameters used for the get_heat_tree_plots function.
#' @param obj A metacoder object.
#' @param title The title used in the heat_tree plot.
#' @param treatment_no The number of treatment groups in the data.
#' @param ... Any of the heat tree parameters list below can be used to change the way the heat_tree
#' output is displayed.  However, this function acts as a default list of parameters.  The memebers
#' of the default list will be overridden by the dot parameters.  Any variable in obj$data$stats_tax_data
#' can be used to manipulate the heat tree parameters.  Function calls from the taxa package must
#' be done explicitely on the metacoder object.
#' @return A list used with do.call and the metacoder::heat_tree function.
#' @pretty_print TRUE
#' @export
#' @family Visualizations
#' @rdname get_heat_tree_parameters
#' @seealso
#'  \code{\link[metacoder]{heat_tree}}
#'
#'  \code{\link[taxa]{n_obs}},  \code{\link[taxa]{taxon_names}}
#'
#'  \code{\link[purrr]{list_modify}}
#'
#'  \code{\link[rlang:quotation]{enquos}},  \code{\link[rlang:quosure]{is_quosure}},  \code{\link[rlang]{eval_tidy}}
#' @importFrom taxa n_obs taxon_names
#' @importFrom purrr list_modify
#' @importFrom rlang enquos is_quosure eval_tidy
get_heat_tree_parameters <- function(obj, title, treatment_no, ...) {
  # Clone the taxmap object.
  input <- obj$clone()

  # Get data used in default params
  # Take advantage of the taxamap internal functions to get the data
  # Use an internal function to get data from the default parameters.
  data_fun <- function(obj, treatment_no, ...) {
    # Get data from default parameters
    data <- obj$data_used(...)
    if (treatment_no == 2) {
      return(list(data = data, params = dplyr::enquos(...)))
    } else if (treatment_no > 2) {
      return(list(data = data, params = rlang::enexprs(...)))
    }
  }

  if (treatment_no > 2) {
    default_parameters <-  data_fun(treatment_no = treatment_no,
      obj = input,
      title = title,
      # NODE
      ## The node size is relevant to the Abundance level
      ## The node color is relevant to wheather the abundance is higher in treatment_1 vs treatment_2 animals
      ## The node labels are relevant to significant taxon names.
      node_size = n_obs,
      node_color = log2_mean_ratio,
      node_label = ifelse(wilcox_p_value < 0.05, taxon_names, NA),
      node_label_size = 1,
      ### The color red indicates higher abundance in Treatment_1 animals
      ### The color blue indicates higher abundance in Treatment_2 animals
      ### The color grey represents no difference in treatment_1 vs treatment_2 abundance
      ### Colorblind node_color_range = c("navy", "grey80", "greenyellow"),
      ### Colorblind node_color_range = c("navy", "grey80", "yellow3"),
      # #ffffbf
      node_color_range = c("#3288bd", "#f1f1f1", "#d53e4f"),
      node_color_trans = "linear",
      node_color_interval = c(-4, 4),
      node_size_axis_label = "Number of OTUs",
      node_color_axis_label = glue::glue("Upregulated (red) vs.\n Downregulated (blue)\n"),
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
      edge_size_axis_label = "Number of OTUs",
      edge_color_axis_label = "Significant Changes \nAmong Treatments",
      # PLOT Options
      initial_layout = "reingold-tilford",
      layout = "davidson-harel",
      repel_labels = TRUE,
      repel_force = 3,
      overlap_avoidance = 3,
      make_edge_legend = FALSE
    )
    #param_list[["data"]] <- "statistical_data"
    param_list <- default_parameters

    # Take advantage of the taxamap internal functions to get the data
    data_list <- input$data_used(...)
    data_list <- c(data_list, default_parameters$data)
    data_list <- data_list[!duplicated(data_list)]

    # Do some cringe worthy magic to dynamically parse parameters
    params <- rlang::enexprs(...) # THese are quosures
    def_param <- list(params = default_parameters$params) # These are quosures
    # Merge lists and replace any matches to the ... args
    param_list <- purrr::list_modify(def_param, params = params)
    param_list <- param_list$params

  } else if (treatment_no == 2 ) {
    default_parameters <-  data_fun(
      obj = obj,
      treatment_no = treatment_no,
      .input = input,
      title = title,
      # NODE
      ## The node size is relevant to the Abundance level
      ## The node color is relevant to wheather the abundance is higher in treatment_1 vs treatment_2
      ## The node labels are relevant to significant taxon names.
      node_size = n_obs,
      node_color = log2_mean_ratio,
      node_label = ifelse(wilcox_p_value < 0.05, taxon_names, NA),
      node_label_size = 1,
      ### The color red indicates higher abundance in Treatment_1 animals
      ### The color blue indicates higher abundance in Treatment_2 animals
      ### The color grey represents no difference in Treatment_1 vs Treatment_2 abundance
      ### Colorblind node_color_range = c("navy", "grey80", "greenyellow"),
      ### Colorblind node_color_range = c("navy", "grey80", "yellow3"),
      # #ffffbf
      node_color_range = c("#3288bd", "#f1f1f1", "#d53e4f"),
      node_color_trans = "linear",
      node_color_interval = c(-4, 4),
      node_size_axis_label = "Number of OTUs",
      node_color_axis_label = "Upregulated (red) vs.\n Downregulated (blue)\n",
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
      edge_size_axis_label = "Number of OTUs",
      edge_color_axis_label = "Significant Changes \nAmong Treatments",
      # PLOT Options
      initial_layout = "reingold-tilford",
      layout = "davidson-harel",
      repel_labels = TRUE,
      repel_force = 3,
      overlap_avoidance = 3,
      make_edge_legend = FALSE
    )
    # Take advantage of the taxamap internal functions to get the data
    data_list <- input$data_used(...)
    data_list <- c(data_list, default_parameters$data)
    data_list <- data_list[!duplicated(data_list)]

    # Do some cringe worthy magic to dynamically parse parameters
    params <- dplyr::enquos(...) # THese are quosures
    def_param <- list(params = default_parameters$params) # These are quosures
    # Merge lists and replace any matches to the ... args
    param_list <- purrr::list_modify(def_param, params = params)
    param_list <- param_list$params

    # Find the quosures and evaluate them
    for (item in names(param_list)) {
      if (rlang::is_quosure(param_list[[item]])) {
        # Use the data_list to evaluae the parameter list.
        param_list[[item]] <- rlang::eval_tidy(param_list[[item]], data = data_list)
      }
    }
  }
  return(param_list)
}


#' @title Save Heat Tree Plots
#' @description This function saves heat tree plots storred in a list object to an output folder.
#' @param htrees A named list of heat trees.
#' @param format The format of the output image.  Default: 'tiff'
#' @param start_path The starting path of the output directory.  Default: 'output'
#' @param ... An optional list of parameters to use in the output_dir function.
#' @return An output directory that contains heat tree plots.
#' @pretty_print TRUE
#' @details This function creates an appropriate output directory, where it saves publication ready
#' plots.
#' @examples
#' \dontrun{
#' if(interactive()){
#' # This example uses data that are no longer available in the MicrobiomeR package,
#' # however, they can be easily generated with \code{\link{MicrobiomeR}{as_analyzed_format}}.
#' library(MicrobiomeR)
#' analyzed_silva <- as_MicrobiomeR_format(MicrobiomeR::raw_silva_2, "analyzed_format")
#' h_trees <- get_heat_tree_plots(analyzed_silva, rank_list = c("Phylum", "Class"))
#' # Save to \emph{./output/heat_trees} folder.
#' save_heat_tree_plots(h_trees)
#'  }
#' }
#' @export
#' @family Visualizations
#' @rdname save_heat_tree_plots
#' @seealso
#'  \code{\link[MicrobiomeR]{output_dir}}
#'
#'  \code{\link[ggplot2]{ggsave}}
#' @importFrom ggplot2 ggsave
#' @importFrom crayon yellow green
#' @importFrom glue glue
save_heat_tree_plots <- function (htrees, format="tiff", start_path = "output", ...) {
  # Create the relative path to the heat_tree plots.  By default the path will be <pwd>/output/<experiment>/heat_trees/<format(Sys.time(), "%Y-%m-%d_%s")>
  # With the parameters set the full path will be <pwd>/output/<experiment>/heat_trees/<extra_path>.
  full_path <- output_dir(start_path = start_path, plot_type = "heat_trees", ...)
  message(glue::glue(crayon::yellow("Saving Heat Trees to the following directory: \n", "\r\t{full_path}")))
  # Iterate the heat_tree plot list and save them in the proper directory
  for (rank in names(htrees)) {
    if (rank != "metacoder_object"){
      message(crayon::green("Saving the {rank} Heat Tree."))
      ggplot2::ggsave(filename = sprintf("%s.heat_tree.%s", rank, format), plot = htrees[[rank]], device = format, path = full_path, dpi = 500)
    }
  }
}
