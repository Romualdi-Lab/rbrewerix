#' Sub sample BrewLOI object
#'
#' Selection can be made on samples or genes or both
#'
#' @param brew a BrewLOI object
#' @param sub_samples character vector with the samples to keep
#' @param sub_genes character vector with the genes to keep
#'
#' @return a BrewLOI object
#' @export
#'
subset <- function(brew, sub_samples=NULL, sub_genes=NULL) {
  checkmate::assertClass(brew, "BrewLOI")

  snp_count <- brew@snp_count
  bi_score <- brew@bi_score
  samples <- brew@samples
  genes <- brew@genes

  if (!is.null(sub_samples)) {
    selection <- samples %in% sub_samples
    if (!sum(selection))
      stop("sub_samples intersection is empty")
    snp_count <- snp_count[, selection]
    bi_score <- bi_score[, selection]
    samples <- samples[selection]
  }

  if (!is.null(sub_genes)) {
    selection <- genes %in% sub_genes
    if (!sum(selection))
      stop("sub_genes intersection is empty")
    snp_count <- snp_count[selection,]
    bi_score <- bi_score[selection,]
    genes <- genes[selection]
  }

  new("BrewLOI", snp_count = snp_count,
      bi_score = bi_score,
      samples = samples,
      genes = genes)
}
