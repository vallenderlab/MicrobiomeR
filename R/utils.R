#' @title Mock Excel "VlookUp" Function
#' @description A function that mimicks excels vlookup, but for R's dataframe.
#' @param lookup_data A vector of items to look up.
#' @param df The dataframe to search.
#' @param match_data The column name to search in the dataframe.
#' @param return_data The column of data to return when matched.
#' @return A vector that contains the items of interest.
#' @pretty_print TRUE
#' @details A function that works like the VLOOKUP function in excel.  This function was
#' borrowed from \url{https://www.r-bloggers.com/an-r-vlookup-not-so-silly-idea/}.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Data Manipulators
#' @rdname vlookup
vlookup <- function(lookup_data, df, match_data, return_data) {
  # TODO: Update the way this returns data.  Allow it to add data to a new column.
  m <- match(lookup_data, df[[match_data]])
  df[[return_data]][m]
}


set_path <- function() {
  library(rstudioapi, quietly = TRUE)
  # Initialize the Environment
  current_path <- getActiveDocumentContext()$path
  # The next line set the working directory to the relevant one:
  setwd(dirname(current_path))
}



#' @title Make a Directory
#' @description This is a custom function for creating a directory.
#' @param dirname The name of the directory to create
#' @param path The path to create the directory in.
#' @return Returns the absolute path of the directory.
#' @pretty_print TRUE
#' @details This function will create a directory only if it doesn't exist
#' and then return the path of that directory.
#' @export
#' @rdname mkdir
mkdir <- function(dirname, path=NULL) {
  # Create the path for the directory
  if (is.null(path)){
    path <- getwd()
  }
  dirpath <- file.path(path, dirname)
  # Create and return the directory path
  if (dir.exists(dirpath)) {
    return(dirpath)
  } else {
    dir.create(dirpath)
    return(dirpath)
  }
}

#' @title Viridis Palette Function
#' @description A function that returns a color palette function based off of the veridis package.
#' @param viridis_number The total number of colors used to generate the viridis palette. Default: 800
#' @param viridis_range Tne range of colors in the viridis palette to use. Default: c(300, 800)
#' @param magma_number The total number of colors used to generate the magma palette. Default: 500
#' @param magma_range The range of colors in the magma palette to use. Default: c(0, 500)
#' @return The output of this function is another function, which takes a number to generate a color
#' palette as a character vector.
#' @details The primary purpose of this function is to return a palette-function for generating virdis style
#' color palettes.  By taking the viridis::viridis() and the viridis::magma() colors, and manipulating
#' them, this function can help create a unique set of colors that you can distinguish on a busy plot.
#' The hopes of this function is to help improve plots that use more than 20 colors.  Use the provided
#' example to view the color palette.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Use the default values
#'  > pal_func <- viridis_palette_func()
#'
#'  # Get a palette with 20 colors
#'  > pal <- pal_func(20)
#'
#'  # Make a pie plot of the colros.
#'  > pie(rep(1, length(pal)), col=pal)
#'  }
#' }
#' @export
#' @family Color Palettes
#' @rdname virdis_palette_func
#' @seealso
#'  \code{\link[viridis]{reexports}}
#'  \code{\link[MicrobiomeR]{get_color_palette}}
#' @importFrom viridis viridis magma
virdis_palette_func <- function(viridis_number=800, viridis_range=c(300,800), magma_number=500, magma_range=c(0, 500)) {
  v_min = viridis_range[1]
  v_max = viridis_range[2]
  m_min = magma_range[1]
  m_max = magma_range[2]
  colorRampPalette(
    unique(c(
      rev(viridis::viridis(viridis_number)[v_min:v_max]), viridis::magma(magma_number)[m_min:m_max])
      )
    )
}

#' @title Get Color palette
#' @description Get a color palette with a specific number of colors.
#' @param pal_func A function that returns the output from grDevoces::colorRampPalette. Default: virdis_palette_func
#' @param color_no The number of colors in the palette. Default: 20
#' @param display Boolean for displaying a pie chart of the palette. Default: TRUE
#' @param ... Parameters for the pal_func.
#' @return Returns a color palette in the form of a character vector.
#' @pretty_print TRUE
#' @details This function is meant to be a plugin style function for user created palettes.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Color Palettes
#' @rdname get_color_palette
get_color_palette <- function(pal_func=virdis_palette_func, color_no=20, display=TRUE, ...) {
  func_params <- list(...)
  pal_func <- pal_func(func_params)
  pal <- pal_func(color_no)
  if (display){
    pie(rep(1, length(pal)), col=pal)
  }
  return(pal)
}


#' @title Object handler
#' @description A function that handles the conversion of objects to metacoder (taxa::taxmap) objects.
#' @param obj An object that contains the data being analyzed.  Can be one of the following:
#' \describe{
#'   \item{Phyloseq Object}{An object generated from the phyloseq package.}
#'   \item{Taxmap Object}{An object generated from the metacoder or taxa package.}
#'   \item{RData file}{An RData file generated from the base::save function.  Can have an extension of .RData or .rda.}
#'   }
#' @return The output generated is a taxmap object.
#' @pretty_print TRUE
#' @details This function is used to convert data to metacoder/taxmap objects for microbiome analysis.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Formatting and Validation
#' @rdname object_handler
#' @seealso
#'  \code{\link[metacoder]{parse_phyloseq}}
#'  \code{\link[tools]{fileutils}}
#' @importFrom metacoder parse_phyloseq
#' @importFrom tools file_ext
object_handler <- function(obj) {
  if (is.null(obj)) {
    stop("Please use a metacoder/phyloseq object or an rdata file.")
  } else {
    if (inherits(obj, "phyloseq")) {
      metacoder_object <- metacoder::parse_phyloseq(obj)
    } else if (inherits(obj, "Taxmap")) {
      metacoder_object <- obj
    } else if (file.exists(obj)) {
      if (tools::file_ext(obj) %in% c("RData", ".rda")) {
        load(file = obj)
        if (!"metacoder_object" %in% ls()) {
          stop("Please provide a loadable .RData/.rda file that contains an object called \"metacoder_object\".")
        }
      }
    }
  }
  return(metacoder_object)
}


#' @title Transposing Tidy Data
#' @description This function transposes tables containing numeric and categorical data using the
#' tidyr package.
#' @param .data A matrix/data_frame/tibble for transposing
#' @param ids The column to transpose by. Default: The first column.
#' @param header_name A name for the numeric data that will be transposed.
#' @param preserved_categories A logical denoting weather categorical data should be conserved.  A
#' value of FALSE will cause all categorical data except the \emph{ids} to be dropped.  A value of
#' TRUE will cause the categorical data to preserved by \emph{tidyr::unite}ing these columns.  Default: TRUE
#' @param separated_categories A vector containing ordered column names to use in a previously transposed
#' and categorically preserved table.  Retransposing with this set should yield an exact replicate of
#' the original data.  Default: NULL
#' @return A transposed data table as a tibble.
#' @pretty_print TRUE
#' @details Transposing can help with preforming operations on the rows of your tibbles.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Data Manipulation
#' @rdname transposer
#' @seealso
#'  \code{\link[tibble]{is_tibble}}
#'  \code{\link[dplyr]{select_all}},\code{\link[dplyr]{select}},\code{\link[dplyr]{reexports}}
#'  \code{\link[tidyr]{gather}},\code{\link[tidyr]{unite}},\code{\link[tidyr]{spread}},\code{\link[tidyr]{separate}}
#'  \code{\link[stringr]{str_detect}},\code{\link[stringr]{str_count}}
#' @importFrom tibble is.tibble
#' @importFrom dplyr select_if select as_tibble
#' @importFrom tidyr gather unite spread separate
#' @importFrom stringr str_detect str_count
transposer <- function(.data, ids = NULL, header_name, preserved_categories = TRUE, separated_categories = NULL) {
  # Verify format
  if (!(is.matrix(.data) | is.data.frame(.data) | tibble::is.tibble(.data))) {
    stop("input not transposable")
  }
  input <- .data
  # Get ids if none are given, defaults to the first column
  if (is.null(ids)) {
    ids <- input[1] %>% colnames()
    warning(sprintf("There were no ids given.  Defaulting to the first column: %s", ids))
  }
  # Get numeric data (columns)
  num_cols <- input %>% dplyr::select_if(is.numeric) %>% dplyr::select_if(!names(.) %in% ids) %>% colnames()

  # Transform
  if (preserved_categories == TRUE) { # All categorical data is preserved
    warning("Categorical data will be united as a string, which can be tidyr::separated after re-transposing.")
    preserved_categories <- input %>% dplyr::select(-one_of(c(num_cols))) %>% colnames()
    trans_data <- input %>%
      dplyr::as_tibble() %>%
      tidyr::gather(key = !!sym(header_name), "_data", -c(preserved_categories)) %>% # Gather columns other that aren't preserved
      tidyr::unite("_categories", c(preserved_categories), sep = "<_>") %>% # Preserve the categorical data in 1 column
      tidyr::spread("_categories", "_data") # Spread categorical data over the numerical data
  } else if (preserved_categories == FALSE) { # Only the ids are preserved
    trans_data <- input %>% select(c(ids, num_cols)) %>%
      dplyr::as_tibble() %>%
      tidyr::gather(key = !!sym(header_name), "_data", -c(ids)) %>%
      tidyr::spread(!!sym(ids), "_data")
  }
  # Look for previously transformed tibbles and separate any united colums
  if (all(stringr::str_detect(trans_data[header_name], "\\<\\_\\>"))) {
    if (!is.null(separated_categories)) {
      trans_data <- trans_data %>% tidyr::separate(col = header_name, into = separated_categories, sep = "<_>")
    } else {
      warning("Separated categories has not been supplied.  Columns will be named as \"category_#\".")
      n_cats <- stringr::str_count(trans_data[[header_name]][1], pattern = "<_>") + 1
      trans_data <- trans_data %>% tidyr::separate(col = header_name, into = paste("category", seq(1:n_cats), sep = "_"), sep = "<_>")
    }
  }
  return(trans_data)
}



#' @title Transforming Tidy Data
#' @description This function transforms a table by rows or by columns.
#' @param .data A matrix/data_frame/tibble for transforming.
#' @param func A function, which can be anonymous, that will be used to transform the data.  (e.g. proportions
#' x/sum(x)).
#' @param by This denotes how the data should be transformed (column/row). Default: 'column'
#' @param ids The column to transpose by. Default: The first column.
#' @param header_name A name for the numeric data that will be transposed.
#' @param preserved_categories A logical denoting weather categorical data should be conserved.  A
#' value of FALSE will cause all categorical data except the \emph{ids} to be dropped.  A value of
#' TRUE will cause the categorical data to preserved by \emph{tidyr::unite}ing these columns.  Default: TRUE
#' @param separated_categories A vector containing ordered column names to use in a previously transposed
#' and categorically preserved table.  Retransposing with this set should yield an exact replicate of
#' the original data.  Default: NULL
#' @param ... Additional arguments passed on to \emph{func}
#' @return A tibble that has been transformed.
#' @pretty_print TRUE
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Data Manipulator
#' @rdname transformer
#' @seealso
#'  \code{\link[MicrobiomeR]{transposer}}
#'  \code{\link[dplyr]{select_all}},\code{\link[dplyr]{select}}
#'  \code{\link[purrr]{modify}}
#' @importFrom dplyr select_if select
#' @importFrom purrr modify_at
transformer <- function(.data, func, by = "column", ids = NULL, header_name = NULL, preserved_categories = TRUE, separated_categories = NULL, ...) {
  if (!(is.matrix(.data) | is.data.frame(.data) | tibble::is.tibble(.data))) {
    stop("input not transposable")
  }
  input <- .data
  if (by == "row") {
    # Transpose table and unite all categorical data into one column for row based transformations
    input <- input %>% MicrobiomeR::transposer(ids = ids, header_name = header_name, preserved_categories = preserved_categories)
  }
  # Get numeric column names
  num_cols <- input %>% dplyr::select_if(is.numeric) %>% dplyr::select_if(!names(.) %in% ids) %>% colnames()
  # Get all other columsn as preserved columns
  preserved_categories <- input %>% dplyr::select(-one_of(num_cols)) %>% colnames()
  # Transform data
  trans_data <- purrr::modify_at(input, num_cols, func, list(...))
  if (by == "row") {
    # Retranspose the table and separate the categorical data for row based transformations
    trans_data <- trans_data %>% MicrobiomeR::transposer(ids = header_name, header_name = ids, separated_categories = separated_categories, preserved_categories = FALSE)
  }
  return(trans_data)
}


#' @title Create an Output Directory
#' @description A function for generating a consistent file system for data files and visualizations.  If none
#' of the parameters are set the full path generated by this function will be \emph{root_path}/output
#' @param start_path  With the start_path set, the full path generated by this function will be
#' \emph{root_path}/\emph{start_path}/\emph{?experiment?}/\emph{?plot_type?}/\emph{end_path}. Default: NULL
#' @param experiment (optional) This will add a second level to the start_path file system. Default: NULL
#' @param plot_type (optional) This will add a third level to the start_path file system.  With \strong{ONLY}
#' the plot_type set the full path generated by this function will be \emph{root_path}/\emph{plot_type}.  Default: NULL
#' @param end_path (optional) This will add the final directory to the start_path file system Default: A formatted date-time string.
#' @param root_path The root of the new file system (start_path or not).  Default: The working directory.
#' @param custom_path This is an absolute path that overrides everything else.  Default: NULL
#' @param overwrite Do you want to overwrite existing output directories?, Default: FALSE
#' @return Creates a directory for output and returns the file.path string.
#' @pretty_print TRUE
#' @details This function is increadibly useful on it's own but also for various other plotting/saving functions within the package.
#'  It helps keep data organized using a standard workflow.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Project Management
#' @rdname get_output_dir
#' @seealso
#'  \code{\link[glue]{glue}}
#'  \code{\link[stringr]{case}}
#' @importFrom glue glue
#' @importFrom stringr str_to_lower
get_output_dir <- function(start_path=NULL, experiment=NULL, plot_type=NULL, end_path=NULL, root_path=NULL, custom_path = NULL, overwrite=FALSE) {
  # Create the relative path to the plots.  By default the full_path will be <root_path>/output
  # With ONLY the plot_type set the full_path will be <root_path>/<plot_type>
  # With the start_path set the full path will be <root_path>/<start_path>/?experiment?/?plot_type?/<end_path>.
  # Experiment and plot_type are optional.  End_path is optional, but defaults to formatted data-time string.
  if (is.null(custom_path)) {
    # Get the root path (generally this already exists; like the working directory)
    if (is.null(root_path)) {
      full_path <- file.path(getwd())
    } else {
      full_path <- file.path(root_path)
    }
    # Add the output folder to use (generally a new folder, or pre-existing output folder)
    if (!is.null(start_path)) {
      full_path <- file.path(full_path, start_path)
      # Add a folder for the specific experiment you are conducting
      if (!is.null(experiment)) {
        full_path <- file.path(full_path, experiment)
      }
      # Add a folder for the plot_type that's being generated
      if (!is.null(plot_type)) {
        full_path <- file.path(full_path, plot_type)
      }
      # Add a folder for any extra
      if (!is.null(end_path)) {
        full_path <- file.path(full_path, extra_path)
      } else if (is.null(end_path)) {
        full_path <- file.path(full_path, format(Sys.time(), "%Y-%m-%d_%s"))
      }
    } else if (is.null(end_path)) {
      # Add a folder for the plot_type that's being generated
      if (!is.null(plot_type)) {
        full_path <- file.path(full_path, plot_type)
      } else {
        full_path <- file.path(full_path, "output")
      }
    }
  } else {
    full_path <- custom_path
  }
  # Directory creation.
  answer_flag <- FALSE
  if (dir.exists(full_path) && overwrite == FALSE) {
    stop(glue::glue("The directory {full_path} already exists. And you don't want to overwrite the directory."))
  } else if (dir.exists(full_path) && overwrite == TRUE) {
    warning(glue::glue("You have chosen to overwrite the directory: {full_path}."))
    while (answer_flag == FALSE) {
      answer <- readline(prompt = "Are you sure? (Y/N)")
      if (stringr::str_to_lower(answer) == "y") {
        dir.create(full_path, recursive = TRUE)
        answer_flag <- TRUE
      } else if (stringr::str_to_lower(answer) == "n") {
        warning(glue::glue("Please set overwrite to FALSE and change the path ({full_path}) with experiment and/or other_path."))
        stop("The files haven't been saved.  You have chosen not to overwrite you files.")
      } else {
        warning("Please enter Y for YES or N for NO.  This is not case sensitive.")
        answer_flag <- FALSE
      }
    }
  } else if (!dir.exists(full_path)) {
    warning(glue::glue("Creating a new directory: {full_path}"))
    dir.create(full_path, recursive = TRUE)
  }
  return(full_path)
}

####################### Private Package Variable #######################
pkg.private <- new.env(parent = emptyenv())

pkg.private$format_level_list <- list(unknown_format = -1, mixed_format = -1,
                                      phyloseq_format = 0, raw_format = 1,
                                      basic_format = 2, analyzed_format = 3)

pkg.private$format_table_list <- list(
  raw_format =  c("otu_abundance", "otu_annotations"),
  basic_format = c("otu_abundance", "otu_annotations", "taxa_abundance", "otu_proportions", "taxa_proportions"),
  analyzed_format = c("otu_abundance", "otu_annotations", "taxa_abundance", "otu_proportions", "taxa_proportions", "statistical_data", "stats_tax_data"),
  phyloseq_format = c("otu_table", "tax_data", "sample_data", "phy_tree"),
  all_formats = c("otu_abundance", "otu_annotations", "taxa_abundance", "otu_proportions", "taxa_proportions", "statistical_data", "stats_tax_data",
                  "otu_table", "tax_data", "sample_data", "phy_tree")
)

pkg.private$ranks <- list(c("Kingdom"),
                          c("Phylum"),
                          c("Class"),
                          c("Order"),
                          c("Family"),
                          c("Genus"),
                          c("Species"))
pkg.private$rank_index <- list(Kingdom = 1,
                               Phylum = 2,
                               Class = 3,
                               Order = 4,
                               Family = 5,
                               Genus = 6,
                               Species = 7)
pkg.private$mc_df_rank_list <- list(Kingdom = 4,
                                    Phylum = 5,
                                    Class = 6,
                                    Order = 7,
                                    Family = 8,
                                    Genus = 9,
                                    Species = 10)
pkg.private$input_files = list(
  biom_files = list(
    silva = "data/silva_OTU.biom",
    greengenes = "data/greengenes_OTU.biom"),
  tree_files = list(
    silva = "data/silva.tre",
    greengenes = "data/greengenes.tre"),
  metadata = "data/nephele_metadata.txt"
)

pkg.private$rdata_files = list(
  phyloseq = list(
    silva = "data/silva_phyloseq_obj.RData",
    greengenes = "data/greengenes_phyloseq_obj.RData"
  ),
  metacoder = list(
    silva = "data/silva_metacoder_obj.RData",
    greengenes = "data/greengenes_metacoder_obj.RData"
  )
)


