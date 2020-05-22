#' For each gene, counts the biallelic samples
#'
#' @param brew a BrewLOI object
#'
#' @export
count_biallelic_samples <- function(brew) {
  checkmate::assertClass(brew, "BrewLOI")
  counts <- apply(brew@bi_score > 0, 1, sum)
  names(counts) <- brew@genes
  counts
}

#' Annot samples
#'
#' @param brew a BrewLOI object
#' @param sample_annot a data frame for sample annotation
#'
#' @export
summarizeBrew <- function(brew, sample_annot) {
  checkmate::assertClass(brew, "BrewLOI")
  bi <- brew@bi_score
  bi[brew@bi_score > 0 ] <- 1
  bi <- t(bi)
  row.names(bi) <- brew@samples
  colnames(bi) <- brew@genes

  if (nrow(sample_annot) != nrow(bi))
    stop("Sample annotation mismatch")

  if (!is.null(row.names(sample_annot)))
    if (!identical(row.names(sample_annot), row.names(bi)))
      stop("sample name mismatch between sample_annotation on brew object")

  cbind(sample_annot, bi)
}
