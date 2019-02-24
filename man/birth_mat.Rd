% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/birth_mat.R
\name{birth_mat}
\alias{birth_mat}
\title{Calculate Births for Each Element of Place of Birth - Place of Residence Stock Matrix}
\usage{
birth_mat(b_por = NULL, m2 = NULL, non_negative = TRUE)
}
\arguments{
\item{b_por}{Vector of numberic values for births in each place of residence}

\item{m2}{Matrix of migrant stock totals at time \emph{t}+1. Rows in the matrix correspond to place of birth and columns to place of residence at time \emph{t}+1.}

\item{non_negative}{Adjust birth matrix calculation to ensure all deductions from \code{m2} will result in positive population counts. On rare occasions when working with international stock data the number of births can exceed the increase in the number of native born population.}
}
\value{
Matrix of place of birth by place of residence for new-born’s
}
\description{
This function is predominantly intended to be used within the ffs routines in the migest package.
}
\seealso{
\code{\link{ipf3_qi}}, \code{\link{ffs_diff}}
}
\author{
Guy J. Abel
}