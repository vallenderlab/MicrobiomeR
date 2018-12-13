#' @title Get Color palette
#' @description Get a color palette with a specific number of colors.
#' @param pal_func A function that returns the output from grDevoces::colorRampPalette. Default: virmag_palette_func
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
get_color_palette <- function(pal_func=virmag_palette_func, color_no=20, display=TRUE, ...) {
  pal_func <- pal_func(...)
  pal <- pal_func(color_no)
  if (display){
    pie(rep(1, length(pal)), col=pal)
  }
  return(pal)
}

#' @title Scico Palette Function
#' @description A function that returns a color palette function based off of the scico package.
#' @param scico_palette The scico palette to use.  Default: 'batlow'
#' @param scico_number The number of colors to use in the scico palette.  Default: 800
#' @param scico_range The range in the color palette to use.  Default: c(0, scico_number)
#' @return The output of this function is another function, which takes a number to generate a color
#' palette as a character vector.
#' @pretty_print TRUE
#' @details The purpose of this function is to provide an interpolated scico palette for using with the
#' get_color_palette function.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Color Palettes
#' @rdname scico_palette_func
#' @seealso
#'  \code{\link[grDevices]{colorRamp}}
#'  \code{\link[scico]{scico}}
#' @importFrom grDevices colorRampPalette
#' @importFrom scico scico
scico_palette_func <- function(scico_palette="batlow", scico_number=800, scico_range=c(0,scico_number)) {
  s_min <- scico_range[1]
  s_max <- scico_range[2]
  crp <- grDevices::colorRampPalette(unique(scico::scico(scico_number, palette = scico_palette)[s_min:s_max]))
  return(crp)
}

#' @title Viridis Palette Function
#' @description A function that returns a color palette function based off of the viridis package.
#' @param viridis_palette The viridis palette to use.  Default: 'batlow'
#' @param viridis_number The number of colors to use in the viridis palette.  Default: 800
#' @param viridis_range The range in the color palette to use.  Default: c(0, viridis_number)
#' @return The output of this function is another function, which takes a number to generate a color
#' palette as a character vector.
#' @pretty_print TRUE
#' @details The purpose of this function is to provide an interpolated viridis palette for using with the
#' get_color_palette function.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Color Palettes
#' @rdname viridis_palette_func
#' @seealso
#'  \code{\link[grDevices]{colorRamp}}
#'  \code{\link[viridis]{viridis}}
#' @importFrom grDevices colorRampPalette
#' @importFrom viridis viridis
viridis_palette_func <- function(viridis_palette="viridis", viridis_number=800, viridis_range=c(0,viridis_number)) {
  s_min <- viridis_range[1]
  s_max <- viridis_range[2]
  crp <- grDevices::colorRampPalette(unique(viridis::viridis(viridis_number, option = viridis_palette)[s_min:s_max]))
  return(crp)
}
combo_palette_func <- function(...) {
  params <- list(...)
  for (item in list) {
    if (inherits(x = item, what = "list")) {

    } else {
      stop("Each parameter you provide must be a \"list\" with 3 variables: \n
           palette:  a character vector of colors (color palette).\n
           number:  an integer used within the ")
    }
  }
}
#' @title Viridis-Magma Palette Function
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
#' @rdname virmag_palette_func
#' @seealso
#'  \code{\link[viridis]{reexports}}
#'  \code{\link[MicrobiomeR]{get_color_palette}}
#' @importFrom viridis viridis magma
virmag_palette_func <- function(viridis_number=800, viridis_range=c(300,800), magma_number=500, magma_range=c(0, 500)) {
  v_min = viridis_range[1]
  v_max = viridis_range[2]
  m_min = magma_range[1]
  m_max = magma_range[2]
  crp <- colorRampPalette(
    unique(
      c(
        rev(viridis::viridis(viridis_number)[v_min:v_max]),
        viridis::magma(magma_number)[m_min:m_max]
      )
    )
  )
  return(crp)
}

