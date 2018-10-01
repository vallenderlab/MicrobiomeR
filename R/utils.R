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
