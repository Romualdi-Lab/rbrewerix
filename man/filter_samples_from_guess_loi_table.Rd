% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_guess_loi_table.R
\name{filter_samples_from_guess_loi_table}
\alias{filter_samples_from_guess_loi_table}
\title{Filter samples of guess loi table using colnames}
\usage{
filter_samples_from_guess_loi_table(
  LOI,
  keep_samples = c(),
  invert_selection = FALSE
)
}
\arguments{
\item{LOI}{a object obtained by read_guess_loi_tavle_v3}

\item{keep_samples}{vector with samples type to keep}

\item{invert_selection}{invert the selection}
}
\value{
a GuessLoi object
}
\description{
Filter samples of guess loi table using colnames
}
