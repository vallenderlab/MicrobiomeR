#' @title Create a table from a data frame
#'
#' @description Creates an image of a table from a dataframe.
#'
#' @param dataframe A dataframe of data.
#'
#' @return The text table as an object that can be saved to png/tiff/jpg/pdf.
#' @importFrom ggpubr ggtexttable ttheme colnames_style tbody_style
#' @export
create_pub_table <- function(dataframe) {
  ggpubr::ggtexttable(dataframe,
    theme = ggpubr::ttheme(
      colnames.style = ggpubr::colnames_style(color = "black", fill = "white", linewidth = 0, linecolor = "white", face = "bold"),
      tbody.style = ggpubr::tbody_style(
        color = "black", face = "plain", size = 12,
        fill = "white", linewidth = 0, linecolor = "white"
      )
    )
  )
}

#' @title Make a Directory
#'
#' @description This is a custom function for creating a directory.
#'
#' @param dirname The name of the directory to create
#' @param path The path to create the directory in.
#'
#' @return Returns the absolute path of the directory.
#' @pretty_print TRUE
#' @details This function will create a directory only if it doesn't exist and then return the path of that directory.
#' @export
#' @rdname mkdir
mkdir <- function(dirname, path = NULL) {
  # Create the path for the directory
  if (is.null(path)) {
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
