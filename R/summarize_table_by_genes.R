#' Summarize Guess LOI Table object by gene
#'
#' @param LOI a object obtained by read_guess_loi_tavle_v3
#' @param thr threshold
#' @param odepth overall depth
#' @param min_ac minimum allele count
#' @param min_number_of_sig_SNP minimal number of bi allelic SNV to call a gene
#'
#' @return a GuessLoi object
#'
#' @export
summarize_table_by_genes <- function(LOI, thr=0.2, odepth=20, min_ac=3, min_number_of_sig_SNP=1) {
  suppressPackageStartupMessages(library(dplyr))
  snvs <- LOI$snv_gene
  ar <- LOI$ar
  pvalues <- LOI$pvalues
  ref_alt_counts <- LOI$ref_alt_counts

  ref <- apply(ref_alt_counts, 2, get_ref_counts)
  alt <- apply(ref_alt_counts, 2, get_alt_counts)
  depth <- ref+alt

  is_dept_enougth <- na2false(depth >= odepth)
  is_biallelic <- na2false(ref >= min_ac & alt >= min_ac & ar >= thr & is_dept_enougth)
  ar[!is_biallelic] <- 0
  ar[!is_dept_enougth] <- NA

  ar_for_aar <- ar
  ar_for_aar[!is_biallelic] <- NA

  snv_gene_ar <- data.frame(symbol=snvs$symbol, ar_for_aar, check.names = F)
  snv_gene_is_bi <- data.frame(symbol=snvs$symbol, is_biallelic, check.names = F)
  snv_gene_found <- data.frame(symbol=snvs$symbol, is_dept_enougth, check.names = F)

  gene_aar = as.data.frame(snv_gene_ar %>% group_by(symbol) %>% summarize_all(mean, na.rm=T))
  biallelic_snvs <- as.data.frame(snv_gene_is_bi %>% group_by(symbol) %>% summarize_all(sum, na.rm=T))
  all_snvs <- as.data.frame(snv_gene_found %>% group_by(symbol) %>% summarize_all(sum, na.rm=T))
  gene_annotation <- as.data.frame(snvs %>% group_by(symbol) %>% summarize_all(head, 1) %>% select(symbol,chr, TSS, type, source))

  invalid <- biallelic_snvs[,-1] < min_number_of_sig_SNP
  gene_aar_value <- gene_aar[,-1]
  gene_aar_value[invalid] <- 0
  gene_aar <- data.frame(symbol=gene_aar$symbol, gene_aar_value, check.names = F)

  biallelic_snvs_value <- biallelic_snvs[,-1]
  biallelic_snvs_value[invalid] <- 0
  biallelic_snvs <- data.frame(symbol=biallelic_snvs$symbol, biallelic_snvs_value, check.names = F)

  if (all (gene_aar$symbol == biallelic_snvs$symbol &
           biallelic_snvs$symbol == all_snvs$symbol &
           all_snvs$symbol == gene_annotation$symbol )) {
    gene_aar$symbol <- NULL
    biallelic_snvs$symbol <- NULL
    all_snvs$symbol <- NULL
  }

  list(annotation=gene_annotation,
       gene_aar=gene_aar,
       gene_all_snvs = all_snvs,
       gene_bi_snvs = biallelic_snvs
  )
}


na2zero <- function(x) {
  x[is.na(x)] <- 0
  x
}


#' Summarize Guess LOI Table object by gene
#'
#' @param LOI a object obtained by read_guess_loi_tavle_v3
#' @param thr threshold
#' @param odepth overall depth
#' @param min_ac minimum allele count
#'
#' @return a GuessLoi object
#'
#' @export
summarize_table_by_genes_big <- function(LOI, thr=0.2, odepth=20, min_ac=3) {
  strong_LOI <- remove_weak_snvs(LOI, thr, odepth, min_ac)
  snvs <- strong_LOI$snv_gene
  ar <- strong_LOI$ar
  ref_alt_counts <- strong_LOI$ref_alt_counts

  gene_wise_list <- tapply(seq_len(nrow(snvs)), snvs$symbol, function(idx) {
    gene_table <- snvs[idx, , drop=F]
    gene_ar <- ar[idx, ,drop=F]
    gene_ref_alt_counts <- ref_alt_counts[idx, ,drop=F]

    gene_ref_alt_counts <- apply(gene_ref_alt_counts, 2, function(x) {
      mat <- do.call(rbind, strsplit(x, ","))
      alt <- sapply(mat[,2], function(x) {ifelse(x=="NA", NA, as.numeric(x))})
      ref <- sapply(mat[,1], function(x) {ifelse(x=="NA", NA, as.numeric(x))})
      mat <- cbind(ref,alt)
      paste(colSums(mat, na.rm=T), collapse=",")
    })

    ann <- gene_table[1, c("symbol", "chr", "TSS", "type", "source")]


    biallelic_snv <- apply(gene_ar>=thr, 2, sum, na.rm=T)

    bi_snv = biallelic_snv
    all_snv <- apply(!is.na(gene_ar), 2, sum)
    average_allelic_ratio <- apply(gene_ar, 2, function(ar) {
      aar <- nan2zero(mean(ar[ar>=thr], na.rm = T))
    })
    # If no biallelic SNVs are found assign all SNV profiles
    average_allelic_ratio[biallelic_snv==0] <- 0
    biallelic_snv[biallelic_snv==0] <- all_snv[biallelic_snv==0]
    bi_ar <- paste(biallelic_snv, average_allelic_ratio, sep=":")
    names(bi_ar) <- names(biallelic_snv)

    list(
      ann=ann,
      aar=data.frame(t(average_allelic_ratio), check.names = F),
      bi_ar=data.frame(t(bi_ar), check.names = F),
      evaluated_snv=data.frame(t(biallelic_snv), check.names = F),
      all_snv=data.frame(t(all_snv), check.names = F),
      bi_snv=data.frame(t(bi_snv), check.names = F),
      gene_ref_alt_counts=data.frame(t(gene_ref_alt_counts), check.names = F)
      )
  })

  gene_annotation <- do.call(rbind, lapply(gene_wise_list, get_gene_annotation))
  row.names(gene_annotation) <- NULL
  gene_annotation$TSS <- as.numeric(as.character(gene_annotation$TSS))

  gene_aar <- do.call(rbind, lapply(gene_wise_list, get_gene_aar))
  row.names(gene_aar) <- NULL

  gene_evaluated_snv <- do.call(rbind, lapply(gene_wise_list, get_gene_e_snv))
  row.names(gene_evaluated_snv) <- NULL

  gene_bi_aar <- do.call(rbind, lapply(gene_wise_list, get_gene_bi_aar))
  row.names(gene_bi_aar) <- NULL

  gene_ref_alt_counts <- do.call(rbind, lapply(gene_wise_list, get_gene_ref_alt_counts))
  row.names(gene_ref_alt_counts) <- NULL

  ## Consider adding tydyr::gather
  ## tidyr::gather(summarized_tables$gene_bi_ar,
  ## "sample", "snvs_ar", -c("symbol", "chr", "TSS", "type", "source"))

  # list(gene_bi_ar=gene_wise_bi_ar, gene_ar=gene_wise_ar, gene_bi=gene_wise_bi)
  list(annotation=gene_annotation,
       gene_aar=gene_aar,
       gene_esnv=gene_evaluated_snv,
       gene_bi_aar = gene_bi_aar,
       gene_ref_alt_counts= gene_ref_alt_counts
       # ,gene_wise_list = gene_wise_list)
  )
}

nan2zero <- function(x) {
  if (is.na(x))
    return(0)
  return(x)
}

get_gene_annotation <- function(gene_obj) {
  gene_obj$ann
}

get_gene_aar <- function(gene_obj) {
  gene_obj$aar
}

get_gene_e_snv <- function(gene_obj) {
  gene_obj$evaluated_snv
}

get_gene_bi_aar <- function(gene_obj) {
  gene_obj$bi_ar
}

get_gene_ref_alt_counts <- function(gene_obj) {
  gene_obj$gene_ref_alt_counts
}
