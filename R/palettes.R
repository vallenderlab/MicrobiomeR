#' @title Get a Color Palette
#' @description Get a color palette with a specific number of colors.
#' @param pal_func A function that \emph{returns the output function} from \code{\link[grDevices]{colorRamp}} or
#' a \emph{the output of the output function} from \code{\link[grDevices]{colorRamp}}.  Default: viridis_magma_palette
#' @param color_no The number of colors in the palette.  Default: 20
#' @param display A logical for displaying a pie chart of the palette. Default: TRUE
#' @param ... Parameters for the \emph{pal_func}.
#' @return Returns a color palette in the form of a character vector.
#' @pretty_print TRUE
#' @details This function is meant to be a plugin style function for users to created palettes.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  library(MicrobiomeR)
#'  data <- MicrobiomeR::analyzed_silva$data$stats_tax_data$Phylum
#'  data_len <- length(unique(data))
#'  custom_pal <- get_color_palette(color_no=data_len)
#'  print(custom_pal)
#'  }
#' }
#' @export
#' @family Color Palettes
#' @rdname get_color_palette
#' @importFrom graphics pie
get_color_palette <- function(pal_func=viridis_magma_palette, color_no=20, display=TRUE, ...) {
  err_handle <- try({
    pal_func <- pal_func(...)
    pal <- pal_func(n=color_no, ...)
  }, silent = TRUE)
  if (inherits(err_handle, "try-error")) {
    pal <- pal_func(n=color_no, ...)
  }
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
#' @return The output of this function is another function (grDevoces::colorRampPalette), which takes
#' a number to generate an interpolated color palette as a character vector.
#' @pretty_print TRUE
#' @details The purpose of this function is to provide an interpolated scico palette for using with the
#' get_color_palette function.
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(Microbiome)
#' basic_pal <- scico_palette(scico_palette="hawaii")(25)
#' better_pal <- get_color_palette(pal_func = scico_palette, color_no = 25)
#'  }
#' }
#' @export
#' @family Color Palettes
#' @rdname scico_palette
#' @seealso
#'  \code{\link[grDevices]{colorRamp}}
#'  \code{\link[scico]{scico}}
#' @importFrom grDevices colorRampPalette
#' @importFrom scico scico
scico_palette <- function(scico_palette="batlow", scico_number=800, scico_range=0:scico_number) {
  pal <- scico::scico(scico_number, palette = scico_palette)[scico_range]
  pal <- unique(pal)
  crp <- grDevices::colorRampPalette(pal)
  return(crp)
}

#' @title Viridis Palette Function
#' @description A function that returns a color palette function based off of the viridis package.
#' @param viridis_palette The viridis palette to use.  Default: 'viridis'
#' @param viridis_number The number of colors to use in the viridis palette.  Default: 800
#' @param viridis_range The range in the color palette to use.  Default: c(0, viridis_number)
#' @return The output of this function is another function (grDevoces::colorRampPalette), which takes
#' a number to generate an interpolated color palette as a character vector.
#' @pretty_print TRUE
#' @details The purpose of this function is to provide an interpolated viridis palette for using with the
#' get_color_palette function.
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(Microbiome)
#' basic_pal <- viridis(viridis_palette="magma")(25)
#' better_pal <- get_color_palette(pal_func = viridis_palette, color_no = 25)
#'  }
#' }
#' @export
#' @family Color Palettes
#' @rdname viridis_palette
#' @seealso
#'  \code{\link[grDevices]{colorRamp}}
#'  \code{\link[viridis:reexports]{viridis}}
#' @importFrom grDevices colorRampPalette
#' @importFrom viridis viridis
viridis_palette <- function(viridis_palette="viridis", viridis_number=800, viridis_range=0:viridis_number) {
  pal <- viridis::viridis(viridis_number, option = viridis_palette)[viridis_range]
  pal <- unique(pal)
  crp <- grDevices::colorRampPalette()
  return(crp)
}


#' @title Combine Color Palettes
#' @description This function uses dynamic arguments (...) in order to combine multiple
#' color palettes together.
#' @param ... You can use any name for your arguments, but the values must be a named list.
#' The list can only have 4 named members:
#' \describe{
#' \item{palette}{This is a palette function that returns a vector of colors.}
#' \item{args}{This is another named list used for the palette function parameters.}
#' \item{range}{This is a range \emph{(1:10)} used to subset the color palette vector.}
#' \item{rev}{This is a logical \emph{(TRUE/FALSE)} used to reverse the color palette.}
#' }
#' You can add as many parameters you want in order to combine as many color palettes
#' as you want.
#' @return The output of this function is another function (grDevoces::colorRampPalette), which takes
#' a number to generate an interpolated color palette as a character vector.
#' @pretty_print TRUE
#' @details This function allows you to combine a varying number of color palettes and gives you
#' the ability to subset and reverse the palettes that are supplied.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Color Palettes
#' @rdname combination_palette
#' @seealso
#'  \code{\link[grDevices]{colorRamp}}
#' @importFrom grDevices colorRampPalette
combination_palette <- function(...) {
  # Set up parameters and variables for the palette functions
  params <- list(...)
  sub_params <- list()
  pal_build <- c()
  # Loop through parameter list that consists of palette functions to combine
  for (item in names(params)) {
    if (inherits(x = params[[item]], what = "list")) {
      sub_params <- params[[item]]
      # Get the required palette function
      pal_func <- sub_params$palette
      # Get the required palette args
      pal_args <- sub_params$args
      # Call the palette function
      pal <- do.call(what = pal_func, args = pal_args)
      # Optionally subset the color palette
      if ("range" %in% names(sub_params)) {
        pal <- pal[sub_params$range]
      }
      # Optionally reverse the color palette
      if ("rev" %in% names(sub_params) && sub_params[["rev"]] == TRUE) {
        pal <- rev(pal)
      }
    } else {
      stop("Each parameter you provide must be a \"list\" with 3 variables: \n

           palette:  a function that creates a palette as a vector of colors.
           args:  a list used with the palette function.
           range (optional):  a ranged used with '[]' to subset the color palette character vector.\n
           rev (optional):  a logical indicating weather or not to reverse the color palette")
    }
    pal_build <- c(pal_build, pal)
  }
  crp <- grDevices::colorRampPalette(unique(pal_build))
  return(crp)
}


#' @title Viridis-Magma Palette Function
#' @description A function that returns a color palette function based off of the veridis package.
#' @param viridis_number The total number of colors used to generate the viridis palette.  Default: 800
#' @param viridis_range Tne range of colors in the viridis palette to use.  Default: 300:viridis_number
#' @param viridis_rev A logical for reversing the viridis palette.  Default: TRUE
#' @param magma_number The total number of colors used to generate the magma palette.  Default: 500
#' @param magma_range The range of colors in the magma palette to use.  Default: 0:magma_number
#' @param magma_rev A logical for reversing the magma palette.  Default: FALSE
#' @param ... These dots are optionally used as both the magma and viridis function parameters.
#' @return The output of this function is another function (grDevoces::colorRampPalette), which takes
#' a number to generate an interpolated color palette as a character vector.
#' @pretty_print TRUE
#' @details The primary purpose of this function is to return a palette-function for generating virdis style
#' color palettes.  By taking the viridis::viridis() and the viridis::magma() colors, and manipulating
#' them, this function can help create a unique set of colors that you can distinguish on a busy plot.
#' The hopes of this function is to help improve plots that use more than 20 colors.  Use the provided
#' example to view the color palette.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Use the default values
#'  > pal_func <- viridis_palette()
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
#' @rdname viridis_magma_palette
#' @seealso
#'  \code{\link[MicrobiomeR]{combination_palette}}
#'  \code{\link[viridis:reexports]{viridis}},  \code{\link[viridis:reexports]{magma}}
#'  \code{\link[grDevices]{colorRamp}}
#' @importFrom viridis viridis magma
viridis_magma_palette <- function(viridis_number = 800, viridis_range = 300:viridis_number, viridis_rev = TRUE,
                                magma_number = 500, magma_range = 0:magma_number, magma_rev = FALSE, ...) {
  if (!missing(...)){
    v_args = list(n=viridis_number, ...)
    m_args = list(n=magma_number, ...)
  } else {
    v_args = list(n=viridis_number)
    m_args = list(n=magma_number)
  }
  crp <- combination_palette(viridis =
                              list(palette = viridis::viridis,
                                   args = v_args,
                                   range = viridis_range,
                                   rev = viridis_rev),
                            magma =
                              list(palette = viridis::magma,
                                   args = m_args,
                                   range = magma_range,
                                   rev = magma_rev))
  return(crp)
}

