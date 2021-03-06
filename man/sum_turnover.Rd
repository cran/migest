% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sum_turnover.R
\name{sum_turnover}
\alias{sum_turnover}
\title{Extract total in-, out- and net-migration totals from an origin-destination migration flow matrix or data frame.}
\usage{
sum_turnover(
  m,
  drop_diagonal = TRUE,
  include_net = TRUE,
  orig_col = "orig",
  dest_col = "dest",
  flow_col = "flow",
  type = "internal"
)
}
\arguments{
\item{m}{A \code{matrix} or data frame of origin-destination flows. For \code{matrix} the first and second dimensions correspond to origin and destination respectively. For a data frame ensure the correct column names are passed to \code{orig_col}, \code{dest_col} and \code{flow_col}.}

\item{drop_diagonal}{Logical to indicate dropping of diagonal terms, where the origin and destination are the same, in the calculation of totals. Default \code{TRUE}.}

\item{include_net}{Logical to indicate inclusion of a net migration total column for each region, in addition to the total in- and out-flows. Default \code{TRUE}.}

\item{orig_col}{Character string of the origin column name (when \code{m} is a data frame rather than a \code{matrix})}

\item{dest_col}{Character string of the destination column name (when \code{m} is a data frame rather than a \code{matrix})}

\item{flow_col}{Character string of the flow column name (when \code{m} is a data frame rather than a \code{matrix})}

\item{type}{Character string to indicate if flows are \code{internal} or \code{international} to indicate if to use \code{region}, \code{tot_in}, \code{tot_out} or \code{country}, \code{tot_imm} and \code{tot_emi} in output.}
}
\value{
A \code{tibble} with total in- and out-flows for each region.
}
\description{
Extract total in-, out- and net-migration totals from an origin-destination migration flow matrix or data frame.
}
\examples{
# matrix
dn <- LETTERS[1:4]
m <- matrix(data = c(0, 100, 30, 70, 50, 0, 45, 5, 60, 35, 0, 40, 20, 25, 20, 0),
            nrow = 4, ncol = 4, dimnames = list(orig = dn, dest = dn), byrow = TRUE)
sum_turnover(m)

\dontrun{
# data frame (tidy) foramt
library(tidyverse)

# download Abel and Cohen (2019) estimates
f <- read_csv("https://ndownloader.figshare.com/files/26239945")

# turnover for single period
f \%>\% 
  filter(year0 == 1990) \%>\%
  sum_turnover(flow_col = "da_pb_closed", type = "international")

# turnover for all periods using group_by
f \%>\% 
  group_by(year0) \%>\%
  sum_turnover(flow_col = "da_pb_closed", type = "international")
}   
}
