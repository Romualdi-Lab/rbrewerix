#' From ASER annotated table split per sample information
#'
#' @param data the imported ASER annotated table
#' @param annotations_cols interger; consider annotation from 1 to annotation_cols
#' @param chr chr col
#' @param pos pos col
#' @param rs rs col
#' @param gene gene col
#'
#' @export
#'
gather_info_per_samples <- function(data, annotation_cols=6, chr=1, pos=2, rs=3, gene=6) {

  if (ncol(data) < 4)
    stop("Not enough annotation data")

  if (ncol(data) <= annotation_cols) {
    stop("No samples to process")
  }

  samples <- lapply(colnames(data)[-c(seq_len(annotation_cols))], function(name) {
    sample <- na.omit(data.frame(chr=data[,chr], pos=as.numeric(data[,pos]),
                                 rs = data[,rs], gene=data[,gene], aexp=data[, name], stringsAsFactors = F))
    sample_aexp <- extract_expression(sample$aexp)
    sample_aexp$rs <- sample$rs
    sample_ratio <- compute_ratio(sample_aexp)
    sample_snps <- cbind(sample[,1:4], sample_aexp, sample_ratio)
    sample_snps
  })

  names(samples) <- colnames(data)[-c(1:6)]
  samples
}

extract_expression <- function(str_vect) {
  splitted <- strsplit(str_vect, split=",")
  num_splitted <- lapply(splitted, as.numeric)
  df <- data.frame(do.call(rbind, num_splitted))
  colnames(df) <- c("ref", "alt", "pvalue")
  df
}

compute_ratio <- function(table) {
  # table=comparison
  table <- as.matrix(table[,c(1:2)])
  apply(table, 1, function(x) {
    sorted_x <- sort(x[1:2])
    sorted_x[1]/sorted_x[2]
  })
}

