#' Remove weak SNVs
#'
#' @param LOI a object obtained by read_guess_loi_tavle_v3
#' @param thr threshold
#' @param odepth overall depth
#' @param min_ac minimum allele count
#'
#' @return a GuessLoi object
#'
#' @export
remove_weak_snvs <- function(LOI, thr=0.2, odepth=20, min_ac=3) {
  snvs <- LOI$snv_gene
  ar <- LOI$ar
  pvalues <- LOI$pvalues
  ref_alt_counts <- LOI$ref_alt_counts

  # depth <- apply(samples, 2, get_overall_depth)
  ref <- apply(ref_alt_counts, 2, get_ref_counts)
  alt <- apply(ref_alt_counts, 2, get_alt_counts)
  depth <- ref+alt

  is_dept_enougth <- na2false(depth >= odepth)
  is_biallelic <- na2false(ref >= min_ac & alt >= min_ac & ar >= thr)

  # ref_alt_counts[!biallelic_mask] <- NA
  ar[!is_biallelic] <- 0
  ar[!is_dept_enougth] <- NA
  list(snv_gene=snvs,
       ref_alt_counts=ref_alt_counts,
       ar=ar,
       pvalues=pvalues)
}

get_ref_counts <- function(sample) {
  ref <- gsub("([^,]+),([^,]+)", "\\1", sample, perl=T)
  suppressWarnings(as.numeric(ref))
}

get_alt_counts <- function(sample) {
  alt <- gsub("([^,]+),([^,]+)", "\\2", sample, perl=T)
  suppressWarnings(as.numeric(alt))
}

