% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metacoder-formatting.R
\name{as_MicrobiomeR_format}
\alias{as_MicrobiomeR_format}
\title{As MicrobiomeR Format}
\usage{
as_MicrobiomeR_format(obj, format, ...)
}
\arguments{
\item{obj}{An object to be converted to a Taxmap object with \code{\link[MicrobiomeR]{create_taxmap}}.}

\item{format}{The name of the format to convert the object to.}

\item{...}{An optional list of parameters to use in the as_*_format function specified
by the format parameter.}
}
\value{
A Taxmap object in the specified format.
}
\description{
Converts a Taxmap object to the specified format.
}
\details{
See the \link{MicrobiomeR_Formats} documentation.
}
\seealso{
\code{\link[MicrobiomeR]{which_format}},  \code{\link[MicrobiomeR]{as_raw_format}},  \code{\link[MicrobiomeR]{as_basic_format}},  \code{\link[MicrobiomeR]{as_analyzed_format}},  \code{\link[MicrobiomeR]{as_phyloseq_format}},  \code{\link[MicrobiomeR]{create_taxmap}},  \code{\link[MicrobiomeR]{order_metacoder_data}}

Other Formatting: \code{\link{MicrobiomeR_Formats}},
  \code{\link{as_analyzed_format}},
  \code{\link{as_basic_format}},
  \code{\link{as_custom_format}},
  \code{\link{as_phyloseq_format}},
  \code{\link{as_raw_format}}, \code{\link{melt_taxmap}},
  \code{\link{treatment_matrix}}
}
\concept{Formatting}
