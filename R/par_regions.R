#' PAR regions
#'
#' @param species mmusculus of hsapiens
#'
#' @export
#'
par_regions <- function(species=c("hsapiens", "mmusculus")) {
  species = species[1]

  m38_regions_x <- matrix(c(169969756, 170931299), nrow=1)
  m38_regions_y <- matrix(c(90745845, 91644698), nrow=1)

  hg38_regions_x <- rbind(par1=c(10000,2781479), par2=c(155701383,156030895))
  hg38_regions_y <- rbind(par1=c(10000,2781479), par2=c(56887903,57217415))

  par <- list(hsapiens=list(X=hg38_regions_x, Y=hg38_regions_y),
              mmusculus=list(X=m38_regions_x, Y=m38_regions_y))
  par[[species]]
}
