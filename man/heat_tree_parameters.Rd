% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heat-tree.R
\name{heat_tree_parameters}
\alias{heat_tree_parameters}
\title{Get Heat Tree Parameters}
\usage{
heat_tree_parameters(obj, title, treatment_no, ...)
}
\arguments{
\item{obj}{A Taxmap object.}

\item{title}{The title used in the heat_tree plot.}

\item{treatment_no}{The number of treatment groups in the data.}

\item{...}{Any of the heat tree parameters list below can be used to change the way the heat_tree
output is displayed.  However, this function acts as a default list of parameters.  The members
of the default list will be overridden by the dot parameters.  Any variable in obj$data$stats_tax_data
can be used to manipulate the heat tree parameters.  Function calls from the taxa package must
be done explicitly on the Taxmap object.}
}
\value{
A list used with do.call and the metacoder::heat_tree function.
}
\description{
This function gets the parameters used for the heat_tree_plots function.
}
\seealso{
\code{\link[metacoder]{heat_tree}}

\code{\link[taxa]{n_obs}},  \code{\link[taxa]{taxon_names}}

\code{\link[purrr]{list_modify}}

\code{\link[rlang:quotation]{enquos}},  \code{\link[rlang:quosure]{is_quosure}},  \code{\link[rlang]{eval_tidy}}

Other Visualizations: \code{\link{alpha_diversity_plot}},
  \code{\link{correlation_data}},
  \code{\link{correlation_plots}},
  \code{\link{correlation_plot}},
  \code{\link{heat_tree_plots}},
  \code{\link{ordination_plots}},
  \code{\link{ordination_plot}}, \code{\link{plot_limits}},
  \code{\link{save_alpha_diversity_plots}},
  \code{\link{save_correlation_plots}},
  \code{\link{save_heat_tree_plots}},
  \code{\link{save_ordination_plots}},
  \code{\link{save_stacked_barplots}},
  \code{\link{stacked_barplots}},
  \code{\link{stacked_barplot}},
  \code{\link{top_coefficients_barplot}}
}
\concept{Visualizations}
