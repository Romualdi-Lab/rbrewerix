#' Associate gene to a color
#'   Function specific of human imprinted gene databases
#'
#' @param data_frame_for_dot_plot a data to plot
#' @param gene_annotation annotations
#'
#' @export
#'
color_the_gene_label <- function(data_frame_for_dot_plot, gene_annotation) {
  used_genes <- levels(data_frame_for_dot_plot$gene)
  annotation_for_colors <- gene_annotation[gene_annotation$symbol %in% used_genes, c("symbol", "info")]

  row.names(annotation_for_colors) <- annotation_for_colors$symbol
  annotation_for_colors <- annotation_for_colors[used_genes,]

  geneimprint_db <- grepl("evidence:Imprinted", annotation_for_colors$info)
  santoni_db <- grepl("mode:Santoni", annotation_for_colors$info)
  valentina_db <- grepl("mode:Valentina", annotation_for_colors$info)

  annotation_for_colors$source <- NA
  annotation_for_colors$source[santoni_db] <- "Santoni"
  annotation_for_colors$source[valentina_db] <- "Valentina"
  annotation_for_colors$source[geneimprint_db] <- "Geneimprint"

  annotation_for_colors$source <- factor(annotation_for_colors$source)
  annotation_for_colors
}

#' Create 3 colors
#' @param name set the names of the three colors
#'
#' @export
#'
brew_db_3colors <- function(names=NULL) {
  colors <- rev(RColorBrewer::brewer.pal(3,"Set2"))
  colors[1] <- "black"
  if (!is.null(names))
    names(colors) <- names
  colors
}

#' Associate gene to a color PAR edition
#'   Function specific of human imprinted gene databases
#'
#' @param data_frame_for_dot_plot a data to plot
#' @param gene_annotation annotations
#'
#' @export
#'
color_type_gene_label <- function(data_frame_for_dot_plot, gene_annotation) {
  used_genes <- levels(data_frame_for_dot_plot$gene)
  annotation_for_colors <- gene_annotation[gene_annotation$symbol %in% used_genes, c("symbol", "type")]

  row.names(annotation_for_colors) <- annotation_for_colors$symbol
  annotation_for_colors <- annotation_for_colors[used_genes,]

  par <- grepl("PAR", annotation_for_colors$type)
  pseudo <- grepl("hasPseudoGene", annotation_for_colors$type)

  annotation_for_colors$source <- ""
  annotation_for_colors$source[pseudo] <- "hasPseudoGene"
  annotation_for_colors$source[par] <- "PAR"

  annotation_for_colors$source <- factor(annotation_for_colors$source)
  annotation_for_colors
}
