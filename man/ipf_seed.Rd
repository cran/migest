% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ipf_seed.R
\name{ipf_seed}
\alias{ipf_seed}
\title{Quickly create IPF seed}
\usage{
ipf_seed(m = NULL, R = NULL, n_dim = NULL, dn = NULL)
}
\arguments{
\item{m}{Matrix, Array or NULL to build seed. If NULL seed will be 1 for all elements.}

\item{R}{Number of rows, columns and possibly n_dimensions for seed matrix or array.}

\item{n_dim}{Numeric integer for the number of n_dimensions - 2 for matrix, 3 or more for an array}

\item{dn}{Vector of character strings for n_dimension names}
}
\value{
An \code{array} or \code{matrix}
}
\description{
This function is predominantly intended to be used within the ipf routines in the migest package.
}
\author{
Guy J. Abel
}
