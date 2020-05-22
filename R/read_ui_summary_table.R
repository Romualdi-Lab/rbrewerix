#' Import Guess LOI Table from file
#'
#' @param file the path to the summary table file
#'
#' @return a BrewLoi object
#'
#' @export
read_ui_summary_table <- function(file) {

  data <- read.table(file, header=T, row.names=1, check.names = F,
                     sep="\t", stringsAsFactors = F, quote="\"")

  split_profiles <- lapply(data, function(profile) {
    t(sapply(strsplit(profile, ":"), as.numeric))
  })

  SNP_count_df <- lapply(split_profiles, function(x) {
    x[,1]
  })
  SNP_count_df <- data.frame(SNP_count_df)

  SNP_heterozigosity <- lapply(split_profiles, function(x) {
    x[,2]
  })
  SNP_heterozigosity <- data.frame(SNP_heterozigosity)

  new("BrewLOI",
      snp_count = SNP_count_df,
      bi_score = SNP_heterozigosity,
      samples = colnames(data),
      genes = row.names(data))
}
