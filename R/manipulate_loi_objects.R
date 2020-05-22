#' Re-name samples of guess loi table
#'
#' @param LOI a object obtained by read_guess_loi_tavle_v3
#' @param newNames vector with new samples names
#'
#' @return a GuessLoi object
#'
#' @export
rename_samples <- function(LOI, newNames) {
  snvs <- LOI$snv_gene
  ar <- LOI$ar
  pvalues <- LOI$pvalues
  ref_alt_counts <- LOI$ref_alt_counts

  if (NCOL(ar) != length(newNames)) {
    stop("Names length mismatch. ", NCOL(ar)," samples available, ", length(newNames), " names found")
  }

  if (any(duplicated(newNames)))
    warning("Duplicated names found")

  colnames(ar) <- newNames
  colnames(pvalues) <- newNames
  colnames(ref_alt_counts) <- newNames

  list(snv_gene=snvs,
       ref_alt_counts=ref_alt_counts,
       ar=ar,
       pvalues=pvalues)
}


#' Sort samples of guess loi table
#'
#' @param LOI a object obtained by read_guess_loi_tavle_v3
#' @param sorter vector with sorted samples names or sorted indexes
#'
#' @return a GuessLoi object
#'
#' @export
sort_samples <- function(LOI, sorter) {
  snvs <- LOI$snv_gene
  ar <- LOI$ar
  pvalues <- LOI$pvalues
  ref_alt_counts <- LOI$ref_alt_counts

  if (NCOL(ar) < length(sorter)) {
    stop("Names length mismatch. ", NCOL(ar)," samples available, ", length(sorter), " names found")
  }

  if (any(duplicated(sorter)))
    warning("Duplicated names found")

  if (is.numeric(sorter)) {
    missing_columns <- setdiff(sorter, seq_along(colnames(ar)))
  } else {
    missing_columns <- setdiff(sorter, colnames(ar))
  }

  if (length(missing_columns))
    stop("Some of your specified columns are missing: ", paste(missing_columns, collpase=", "))
  ar = ar[, sorter, drop=F]
  pvalues = pvalues[, sorter, drop=F]
  ref_alt_counts = ref_alt_counts[, sorter, drop=F]

  list(snv_gene=snvs,
       ref_alt_counts=ref_alt_counts,
       ar=ar,
       pvalues=pvalues)
}
