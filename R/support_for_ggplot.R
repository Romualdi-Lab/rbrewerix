#' Remove axis title and legend
#'
#' @return list with ggplot2 specification
#'
#' @export
#'
empty_labs <- function() {
  list(
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none")
  )
}

