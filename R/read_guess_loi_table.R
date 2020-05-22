#' Import Guess LOI Table from file
#'
#' @param file the path to the file.
#'   File must have 9 annotation columns
#'
#' @return a BrewLoi object
#'
#' @export
read_guess_loi_table <- function(file, version=c(4,3,2)) {
  readable_versions = c(4,3,2)
  data <- read.table(file, header=T, check.names = F, sep="\t", stringsAsFactors = F, quote="\"", comment.char = '#')

  version=version[1]
  if (! version %in% readable_versions)
    stop("Unable to read specified version")

  if (version %in% c(4,3) & colnames(data)[9] != "TSS")
    stop("Unexpected table format\ntry read_ui_summary_table function if you file is a summary or version=2")

  if (version == 2 & colnames(data)[7] != "TSS")
    stop("Unexpected table format\ntry read_ui_summary_table function if you file is a summary or version=3")

  annotation_end = 9

  if (version == 2) {
    annotation_end = 7
  }

  annotations <- data[,1:annotation_end]
  samples <- data[,-c(1:annotation_end), drop=F]

  if (version == 2) {
    annotations = cbind(annotations[,1:6,drop=F], type=NA, info=NA, annotations[,7,drop=F])
  }

  pvalue_mat <- apply(samples, 2, get_pvalue)
  allelic_ratio_mat <- apply(samples, 2, get_ratio)
  ref_alt_counts <- apply(samples, 2, get_ref_alt_counts)

  list(snv_gene = annotations, ref_alt_counts=ref_alt_counts, ar=allelic_ratio_mat, pvalues=pvalue_mat)
}
# guess-GSE45719-mouse-sc-easrly-stages


get_pvalue <- function(sample) {
  p <- gsub("([^,]+),([^,]+),([^,]+)", "\\3", sample, perl=T)
  suppressWarnings(as.numeric(p))
}

.get_pvalue <- function(sample) {
  sapply(strsplit(sample, ','), function(x) {
    as.numeric(x[3])
  })
}


.get_ratio <- function(sample) {
  sapply(strsplit(sample, ','), function(x) {
    ac = as.numeric(x[1:2])
    min(ac)/max(ac)
  })
}

get_ratio <- function(sample) {
  ref <- gsub("([^,]+),([^,]+),([^,]+)", "\\1", sample, perl=T)
  alt <- gsub("([^,]+),([^,]+),([^,]+)", "\\2", sample, perl=T)
  ratio <- suppressWarnings(as.numeric(ref)/as.numeric(alt))
  select = ratio > 1 & !is.na(ratio)
  ratio[select] <- ratio[select]^-1
  ratio
}


get_ref_alt_counts <- function(sample) {
  ra = gsub("([^,]+),([^,]+),([^,]+)", "\\1,\\2", sample, perl=T)
  ra[is.na(ra)] <- "NA,NA"
  ra
}

.get_ref_alt_counts <- function(sample) {
  sapply(strsplit(sample, ','), function(x) {
    paste(x[1], x[2], sep=",")
  })
}

get_overall_depth <- function(sample) {
  ref <- gsub("([^,]+),([^,]+),([^,]+)", "\\1", sample, perl=T)
  alt <- gsub("([^,]+),([^,]+),([^,]+)", "\\2", sample, perl=T)
  suppressWarnings(as.numeric(ref)+as.numeric(alt))
}
