#' @title Create a table from a data frame
#'
#' @description Creates an image of a table from a dataframe.
#'
#' @param dataframe A dataframe of data.
#'
#' @return The text table as an object that can be saved to png/tiff/jpg/pdf.
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
