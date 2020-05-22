#' Plor gene summary table
#'
#' @param summary_table a object obtained by summarize_table_by_genes
#'
#' @return NULL
#'
#' @export
#'
plot_gene_average_allelic_ratio <- function(summary_table) {
  require(ggplot2)
  require(tidyr)
  gene_bi_aar <- cbind(summary_table$annotation, summary_table$gene_bi_aar)
  gene_bi_aar <- gene_bi_aar[order(gene_bi_aar$chr, gene_bi_aar$TSS), ]
  gene_order <- gene_bi_aar$symbol

  expanded_summary <- tidyr::gather(gene_bi_aar, "sample", "snvs_aar",
                                    -c("symbol", "chr", "TSS", "type", "info"))
  expanded_values <- do.call(rbind, strsplit(expanded_summary$snvs_aar, ":"))
  snv_number=as.numeric(expanded_values[,1])
  snv_number[snv_number > 6] <- 6
  allelic_ratio=as.numeric(expanded_values[,2])
  expanded_summary$snv_number <- snv_number
  expanded_summary$allelic_ratio <- allelic_ratio
  expanded_summary <- expanded_summary[order(expanded_summary$chr, expanded_summary$TSS), ]
  expanded_summary$symbol <- factor(expanded_summary$symbol, levels=rev(gene_order))

  plot <- ggplot(expanded_summary, aes(x=sample, y=symbol, colour=allelic_ratio)) +
    geom_point(aes(size=snv_number)) +
    scale_color_gradientn(colours = c("#ffffff", "#cd6133", "#cc8e35", "#ffa502","#ffa502"), values = c(0,0.2,0.4,0.8,1)) +
    theme_minimal()
  list(plot=plot, data=expanded_summary)
}
