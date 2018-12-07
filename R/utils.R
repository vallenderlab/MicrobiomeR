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


# A function for transposing tibbles with categorical data.
#' @title Transposer function for tidy data.
#' @description This function transposes data using the tidyr package.
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
    trans_data <- input %>%
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
