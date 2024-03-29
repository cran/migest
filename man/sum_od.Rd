% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sum_od.R
\name{sum_od}
\alias{sum_od}
\title{Extract a classic origin-destination migration flow matrix.}
\usage{
sum_od(x = NULL, zero_diag = TRUE, add_margins = TRUE)
}
\arguments{
\item{x}{Array of origin-destination matrices, where the first and second dimensions correspond to origin and destination respectively. Higher dimension(s) refer to additional migrant characteristic(s).}

\item{zero_diag}{Logical to indicate if to set diagonal terms to zero. Default \code{TRUE}.}

\item{add_margins}{Logical to indicate if to add row and column for immigration and emigration totals. Default \code{TRUE}}
}
\value{
Matrix from summing over the first and second dimension. Set diagonals to zero.

Returns a \code{matrix} object of origin-destination flows
}
\description{
Extract a classic origin-destination migration flow matrix from a more detailed dis-aggregation of flows stored in an (\code{array}).
Primarily intended to work with output from \code{\link{ffs_demo}}.
}
\seealso{
\code{\link{ffs_demo}}
}
