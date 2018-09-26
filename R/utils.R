

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
