% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ASER-table-parser.R
\name{gather_info_per_samples}
\alias{gather_info_per_samples}
\title{From ASER annotated table split per sample information}
\usage{
gather_info_per_samples(
  data,
  annotation_cols = 6,
  chr = 1,
  pos = 2,
  rs = 3,
  gene = 6
)
}
\arguments{
\item{data}{the imported ASER annotated table}

\item{chr}{chr col}

\item{pos}{pos col}

\item{rs}{rs col}

\item{gene}{gene col}

\item{annotations_cols}{interger; consider annotation from 1 to annotation_cols}
}
\description{
From ASER annotated table split per sample information
}
