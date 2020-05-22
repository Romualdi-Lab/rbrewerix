check_BrewLOI <- function(object) {
  if (ncol(object@snp_count) != ncol(object@bi_score)){
    msg <- "snp_count and bi_score matrix must have the same column number"
    return(msg)
  }
  if (nrow(object@snp_count) != nrow(object@bi_score)){
    msg <- "snp_count and bi_score matrix must have the same row number"
    return(msg)
  }

  if (ncol(object@snp_count) != length(object@samples)){
    msg <- "samples length mismatch"
    return(msg)
  }

  if (nrow(object@snp_count) != length(object@genes)){
    msg <- "gene length mismatch"
    return(msg)
  }
  return(TRUE)
}

#' BrewLOI-class
#'
#' Class to store the info in the guess_LOI_table
#'
#' @slot snp_count how many SNPs are detected per gene.
#' @slot bi_score the average biallelic score of the SNPs detected per gene
#' @slot samples the sample names
#' @slot genes the gene names
#'
setClass("BrewLOI", package = "brewLOI",
         slots = c(snp_count = "data.frame",
                   bi_score = "data.frame",
                   samples = "character",
                   genes = "character"),
         validity = check_BrewLOI)
