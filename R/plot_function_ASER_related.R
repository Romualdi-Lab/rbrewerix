#' Plot the Chromosome coverage and False Positives
#'
#' @param asample a samnple derived from gather_info_per_samples
#' @param bin_dim binning dimension to compute the coverage
#' @param ceil_to upper limit for coverage to increase resolution at low coverage
#' @param hetero_thr biallelic threshold
#' @param min_depth min depth
#' @param min_alt_counts minimal alternative counts
#' @param pvalue_thr threshold for pvalue
#' @param chromosomes plot chromosomes
#' @param plot_bars plot bar coverage
#' @param plot_smoothed_cov plot smoothed coverage
#'
#' @export
#'
plot_x_false_positive <- function(asample, bin_dim=1000000, ceil_to = 100,
                                          hetero_thr=0.2, min_depth=20, min_alt_counts=3, pvalue_thr=NULL,
                                          chromosomes=c("X"), plot_bars=FALSE, plot_smoothed_cov=FALSE) {
  library(ggplot2)
  chr = "X"
  x_asample <- asample[asample$chr==chr, ]
  x_asample <- x_asample[order(x_asample$pos),]
  breaks <- seq(1, max(x_asample$pos)+bin_dim+40000, bin_dim)

  x_asample$coverage <- x_asample$ref + x_asample$alt
  clusters <- cut(x_asample$pos, breaks, labels=FALSE)
  # clusters_explained <- cut(x_asample$pos, breaks)

  binned_coverage <- tapply(x_asample$coverage, clusters, mean)
  binned_snp_number <- tapply(x_asample$pos, clusters, length)

  # tapply(x_asample$rs, clusters, function(x)x)[["155"]]

  data_coverage <- data.frame(x=as.numeric(names(binned_coverage)), cov=binned_coverage, type="all")
  data_coverage$cov_ori <- data_coverage$cov
  data_coverage$cov[data_coverage$cov > ceil_to] <- ceil_to

  binned_snp_number <- tapply(x_asample$pos, clusters, length)
  data_snp <- data.frame(x=as.numeric(names(binned_snp_number)), count=binned_snp_number/100)
  # tapply(x_asample$rs, clusters, function(x)x)[["155"]]

  # extract_snps_hetero
  hetero_snps_x <- x_asample[x_asample$sample_ratio > 0,c("pos", "sample_ratio", "coverage")]
  cluster_hetero <- cut(hetero_snps_x$pos, breaks, labels=FALSE)

  binned_het_coverage <- tapply(hetero_snps_x$coverage, cluster_hetero, mean)
  data_het_coverage <- data.frame(x=as.numeric(names(binned_het_coverage)), cov=binned_het_coverage, type="het")

  data_points <- data.frame(x=cluster_hetero, het=as.numeric(hetero_snps_x$sample_ratio))
  data_points$hetero <- data_points$het >= 0.2
  data_points$het <- data_points$het*ceil_to
  sum(data_points$hetero)

  data_coverage$PAR <- FALSE
  data_coverage$PAR[data_coverage$x <= cut(par_regions()$X["par1",2], breaks, labels=FALSE)] <- TRUE
  data_coverage$PAR[data_coverage$x >= cut(par_regions()$X["par2",1], breaks, labels=FALSE)] <- TRUE

  # cut(par_regions()$X["par2",1], breaks, labels=FALSE)

  library(ggplot2)
  p <- ggplot(data_coverage, aes(x,cov))

  if (plot_bars) {
    p <- p + geom_bar(stat="identity", position=position_dodge(), alpha=0.4, color="#596275")
  }

  if (plot_smoothed_cov) {
    p <- p + geom_smooth(aes(x,cov_ori, fill=NULL), se=F, span=0.5, color="#ffda79")
  }

  # scale_fill_manual(values=c("#d2dae2", "#05c46b"))
  # geom_line(data=data_points, aes(x,het)) +
  # stat_smooth(data=data_points, aes(x,het), se=F, span=0.1)

  p <- p + geom_smooth(data=data_snp, aes(x,count, fill=NULL), se=F, span=0.5, color="#ffda79")

  p1 <- p +
    geom_vline(xintercept = ceiling(par_regions()$X["par1",2]/bin_dim)+0.5, color="red", size=1) +
    geom_vline(xintercept = floor(par_regions()$X["par2",1]/bin_dim)+1.5, color="red", size=1) +
    # geom_point(data=data_points, aes(x, het, color=hetero, fill=NULL), alpha=0.5) +
    # scale_color_manual(values=c("#485460", "#ff793f")) +
    theme_classic()
  p1

  thr_df <- extract_snps(asample, hetero_thr, min_depth, min_alt_counts, pvalue_thr=pvalue_thr, chromosomes=chromosomes)
  false_pos <- remove_par_genes(thr_df)
  number_of_snps <- NROW(x_asample[x_asample$coverage > min_depth, ,drop=F])
  number_of_snps_no_par <- NROW(remove_par_genes(x_asample[x_asample$coverage >= min_depth, ,drop=F]))

  data_points_thr<- create_hetero_points(thr_df, breaks)
  data_points_thr$het <- data_points_thr$het*ceil_to

  data_points$hetero <- NULL
  data_points$hetero <- FALSE
  data_points_thr$hetero <- TRUE
  data_points <- rbind(data_points, data_points_thr)

  # p2 <- p1 + geom_point(data=data_points_thr, aes(x, het, fill=NULL), color="#ff793f")
  p2 <- p1 +
    geom_point(data=data_points, aes(x, het, color=hetero, shape=hetero, fill=NULL, alpha=hetero)) +
    scale_shape_manual(values=c(1,19)) +
    scale_alpha_manual(values=c(0.5,0.9))+
    scale_color_manual(values=c("#84817a", "#ff5252"))

  list(plot=p2, fp=NROW(false_pos), all_snp=number_of_snps_no_par, fp_ratio=round(NROW(false_pos)/number_of_snps_no_par*10^6))
}

extract_snps <- function(sample_snps, hetero_thr, min_depth, min_alt_counts, pvalue_thr=NULL, chromosomes=NULL) {
  selected <- sample_snps$sample_ratio >= hetero_thr &
    (sample_snps$alt + sample_snps$ref) >= min_depth &
    (sample_snps$alt >=min_alt_counts & sample_snps$ref >= min_alt_counts)

  if (!is.null(pvalue_thr))
    selected <- selected & sample_snps$pvalue <= pvalue_thr

  if (!is.null(chromosomes))
    selected <- selected & sample_snps$chr %in% chromosomes
  sample_snps[selected, ]
}

remove_par_genes <- function(df, invert=F) {
  par_x1 <- list(chr="X",start=0, end=2781479)
  par_x2 <- list(chr="X",start=155701383, end=156030895)
  par_y1 <- list(chr="Y",start=0, end=2781479)
  par_y2 <- list(chr="Y",start=56887903, end=57217415)
  snp_in_x <- df$chr == "X"
  warning("Par defined by the function non official par used")

  df_pos <- as.numeric(df$pos)

  snps_in_par_x <- (df_pos >= par_x1$start & df_pos <=par_x1$end) |
    (df_pos >= par_x2$start & df_pos <=par_x2$end)

  snp_in_y <- df$chr == "Y"
  df_pos <- as.numeric(df$pos)
  snps_in_par_y <- (df_pos >= par_y1$start & df_pos <=par_y1$end) |
    (df_pos >= par_y2$start & df_pos <=par_y2$end)

  select <- !(snp_in_x & snps_in_par_x) & !(snp_in_y & snps_in_par_y)

  if (invert)
    select <- !select

  df$pos <- as.numeric(df$pos)
  no_par_df <- df[select,  , drop=F]
  no_par_df[order(no_par_df$chr, no_par_df$pos), , drop=F]
}

create_hetero_points <- function(hetero_snps_x, breaks) {
  cluster_hetero <- cut(hetero_snps_x$pos, breaks, labels=FALSE)
  coverage <- hetero_snps_x$ref + hetero_snps_x$alt
  binned_het_coverage <- tapply(coverage, cluster_hetero, mean)
  data_het_coverage <- data.frame(x=as.numeric(names(binned_het_coverage)), cov=binned_het_coverage, type="het")

  data.frame(x=cluster_hetero, het=as.numeric(hetero_snps_x$sample_ratio))
}


