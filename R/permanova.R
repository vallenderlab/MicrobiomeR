#' @title Top Coefficients Barplot
#' @description FUNCTION_DESCRIPTION
#' @param top_coefficients PARAM_DESCRIPTION
#' @return Returns a barplot of the top coefficients
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname topcoef_barplot
topcoef_barplot <- function(top_coefficients) {
  # Set graphical parameters
  par(mar = c(3, 16, 2, 2))

  # Plot the top coefficients
  plot <- barplot(sort(top_coefficients),
    horiz = T, las = 1, main = "Top taxa",
    col = ifelse(sort(top_coefficients) >= 0, "#3288bd", "#d53e4f")
  )
  return(plot)
}
