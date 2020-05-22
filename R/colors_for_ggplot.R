#' Blues palette for cell percentages
#'
#' @param na.value color for na
#'
#' @return scale_fill_gradient object with 6 blues
#'
#' @export
#'
cell_percentage_gradient_blues <- function(na.value = "grey50") {
  colors = RColorBrewer::brewer.pal(9, "Blues")[c(2,3,5:8)]

  greens = RColorBrewer::brewer.pal(9, "Greens")[c(2,3,4)]
  reds = RColorBrewer::brewer.pal(9, "Reds")[c(5:7)]

  list(ggplot2::scale_fill_gradientn(colors=colors, na.value = na.value))
}

#' Green Red palette for cell percentages
#'
#' @param na.value color for na
#'
#' @return scale_fill_gradient object with 3 green bottom and 3 red up
#'
#' @export
#'
cell_percentage_gradient_grReds <- function(na.value = "grey50") {
  greens = RColorBrewer::brewer.pal(9, "Greens")[c(1,3,4)]
  reds = RColorBrewer::brewer.pal(9, "Reds")[c(5:7)]
  colors = c(greens, reds)

  list(ggplot2::scale_fill_gradientn(colors=colors, na.value = na.value))
}
