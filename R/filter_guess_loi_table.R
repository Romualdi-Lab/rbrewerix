#' Summarize Guess LOI Table object by gene
#'
#' @param LOI a object obtained by read_guess_loi_tavle_v3
#' @param keep_type vector with gene type to keep
#' @param keep_sources vector with gene imprinting sources to keep
#' @param invert_selection invert the selection
#'
#' @return a GuessLoi object
#'
#' @export
filter_guess_loi_table <- function(LOI, keep_type=c(), keep_source=c(), invert_selection=FALSE) {
  snvs <- LOI$snv_gene
  ar <- LOI$ar
  pvalues <- LOI$pvalues
  ref_alt_counts <- LOI$ref_alt_counts

  source_filters <- NULL
  if (!is.null(keep_source)) {
    source_filters <- do.call(cbind, lapply(keep_source, function(source) {
      grepl(source, snvs$source)
    }))
  }

  type_filters = NULL
  if (!is.null(keep_type)){
    type_filters <- do.call(cbind,lapply(keep_type, function(type) {
      grepl(type, snvs$type)
    }))
  }

  filters <- cbind(source_filters, type_filters)
  if (is.null(filters) & !invert_selection)
    return(LOI)

  if (is.null(filters) & invert_selection)
    return(NULL)

  keep <- apply(filters, 1, any)
  if (invert_selection)
    keep <- !keep

  list(snv_gene=snvs[keep, ,drop=F],
       ref_alt_counts=ref_alt_counts[keep, , drop=F],
       ar=ar[keep, ,drop=F],
       pvalues=pvalues[keep, , drop=F])
}

#' Filter samples of guess loi table using colnames
#'
#' @param LOI a object obtained by read_guess_loi_tavle_v3
#' @param keep_samples vector with samples type to keep
#' @param invert_selection invert the selection
#'
#' @return a GuessLoi object
#'
#' @export
filter_samples_from_guess_loi_table <- function(LOI, keep_samples=c(), invert_selection=FALSE) {
  snvs <- LOI$snv_gene
  ar <- LOI$ar
  pvalues <- LOI$pvalues
  ref_alt_counts <- LOI$ref_alt_counts

  available_samples <- colnames(ref_alt_counts)

  if (!is.null(keep_samples)){
    keep <- available_samples %in% keep_samples
  }
  if (invert_selection)
    keep <- !keep

  list(snv_gene=snvs,
       ref_alt_counts=ref_alt_counts[,keep , drop=F],
       ar=ar[,keep ,drop=F],
       pvalues=pvalues[,keep , drop=F])
}


#' Extract sexual chromosomes
#'
#' @param LOI a object obtained by read_guess_loi_tavle
#' @param invert_selection revert selection
#' @param chrs vector with chromosomes of interest
#'
#' @return a GuessLoi object
#'
#' @export
extract_chromosomes <- function(LOI, invert_selection = F, chrs=c("X","Y")) {
  snvs <- LOI$snv_gene
  ar <- LOI$ar
  pvalues <- LOI$pvalues
  ref_alt_counts <- LOI$ref_alt_counts

  selection <- snvs$chr %in% chrs
  if (invert_selection)
    selection <- !selection

  if (!sum(selection))
    return(NULL)

  list(snv_gene=snvs[selection, ,drop=F],
       ref_alt_counts=ref_alt_counts[selection, , drop=F],
       ar=ar[selection, ,drop=F],
       pvalues=pvalues[selection, , drop=F])
}


#' Filter tables
#'
#' @param LOI a object obtained by read_guess_loi_tavle_v3
#' @param thr threshold
#' @param odepth overall depth
#' @param min_ac minimum allele count
#'
#' @return a GuessLoi object
.remove_weak_snvs <- function(LOI, thr=0.2, odepth=20, min_ac=3) {
  snvs <- LOI$snv_gene
  ar <- LOI$ar
  pvalues <- LOI$pvalues
  ref_alt_counts <- LOI$ref_alt_counts

  biallelic_mask <- apply(ref_alt_counts, 2, is_biallelic, ar_thr=thr, odepth=odepth, min_ac=min_ac)
  depth_mask <- apply(ref_alt_counts, 2, check_depth_counts, odepth=odepth)

  # ref_alt_counts[!biallelic_mask] <- NA
  ar[!biallelic_mask] <- 0
  ar[!depth_mask] <- NA
  list(snv_gene=snvs,
       ref_alt_counts=ref_alt_counts,
       ar=ar,
       pvalues=pvalues)
}

.is_biallelic <- function(sample_ref_alt_count, ar_thr, odepth, min_ac, sep=",") {
  sapply(strsplit(sample_ref_alt_count, sep), function(ac) {
    if (ac[1]=='NA' & ac[2]=="NA")
      return(FALSE)
    ac_num <- c(as.numeric(ac[1]),as.numeric(ac[2]))
    is_above_odepth(ac_num, odepth) && is_above_min_ac(ac_num, min_ac) && is_above_ar_thr(ac_num, ar_thr)
  })
}

is_biallelic <- function(sample_ref_alt_count, ar_thr, odepth, min_ac) {
  ref <- silent_as.numeric(gsub("([^,]+),([^,]+)", "\\1", sample_ref_alt_count, perl=T))
  alt <- silent_as.numeric(gsub("([^,]+),([^,]+)", "\\2", sample_ref_alt_count, perl=T))
  ratio <- ref/alt
  sum <- ref+alt
  select = ratio > 1 & !is.na(ratio)
  ratio[select] <- ratio[select]^-1
  is_bi <- ratio >= ar_thr & ref >= min_ac & alt >= min_ac & sum >= odepth
  na2false(is_bi)
}

.check_depth_counts <- function(sample_ref_alt_count, odepth, sep=",") {
  sapply(strsplit(sample_ref_alt_count, sep), function(ac) {
    if (ac[1]=='NA' & ac[2]=="NA")
      return(FALSE)
    ac_num <- c(as.numeric(ac[1]),as.numeric(ac[2]))
    is_above_odepth(ac_num, odepth)
  })
}

check_depth_counts <- function(sample_ref_alt_count, odepth) {
  ref <- silent_as.numeric(gsub("([^,]+),([^,]+)", "\\1", sample_ref_alt_count, perl=T))
  alt <- silent_as.numeric(gsub("([^,]+),([^,]+)", "\\2", sample_ref_alt_count, perl=T))
  sum <- ref+alt
  above <- sum >= odepth
  na2false(above)
}

silent_as.numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

na2false <- function(x) {
  x[is.na(x)] <- FALSE
  x
}

is_above_odepth <- function(ac, odepth){
  sum(ac) >= odepth
}

is_above_min_ac <- function(ac, min_ac){
  min(ac) >= min_ac
}

is_above_ar_thr <- function(ac, ar_thr){
  min(ac)/max(ac) >= ar_thr
}
