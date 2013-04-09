\name{ffs}
\alias{ffs}
\title{
Flows from Stocks Methodology
}
\description{
Estimates migrant transitions flows between two sequential migrant stock tables.
}
\usage{
ffs(P1, P2, d, b, m = NULL, net=FALSE, d.adj=FALSE, ...)
}
\arguments{
  \item{P1}{
Matrix of migrant stock totals at time \emph{t}. Rows in the matrix correspond to place of birth and columns to place of residence at time \emph{t} 
}
  \item{P2}{
Matrix of migrant stock totals at time \emph{t}+1. Rows in the matrix correspond to place of birth and columns to place of residence at time \emph{t}+1.
}
  \item{d}{
Number of deaths between time \emph{t} and \emph{t}+1 in each region.
}
  \item{b}{
Number of births between time \emph{t} and \emph{t}+1 in each region.
}
  \item{m}{
Array of auxiliary data. By default set to 1 for all origin-destination-migrant type combinations.
}
  \item{net}{
Calcualate adjusted stock tables to match the row totals of the stock table after demographic accounting. By default \code{FALSE}.
}
  \item{d.adj}{
Calculate deaths by place of residence table to match the row totals of the stock table after demographic accounting. By default \code{FALSE}.
}
  \item{\dots}{
Additional arguments passes to \code{\link{ipf3.qi}}.
}
}
\details{
Estimates migrant transitions flows between two sequential migrant stock tables as shown in Abel (2013). The length of \code{b} and \code{d} must equal the number of rows in \code{P1} and number of columns in \code{P2}.

Setting \code{net=TRUE} uses an alternative demographic accounting method to adjusted stock tables to match the row totals of the stock table after demographic accounting. Setting \code{adj.d=TRUE} uses the calcualtion of the deaths by place of residence table to match the row totals of the stock table after demographic accounting. Both of these options maintain the net migration flow implied by the population, birth and death data. Only one of these arguments can be set to \code{TRUE}. Discussion of these methods is given in Abel (2013b).
}
\value{
Returns a \code{list} object with:
  \item{mu }{Array of indirect estimates of origin-destination matrices by place of birth.}
  \item{it }{Iteration count.}
  \item{tol }{Tolerance level at final iteration.}
  \item{y }{Array of indirect estimates of origin-destination matrices by place of birth with additional rows and columns for births, deaths and moves to other regions.}
}
\references{
Abel, G. J. (2013). Estimating Global Migration Flow Tables Using Place of Birth. \emph{Demographic Research} 28, (18) 505-546

Abel, G. J. (2013b). Combining Bilateral Migrant Stock and Net Migration Flow Information to Estimate Global Flow Tables. \emph{Vienna Institute of Demography Working Paper} Forthcoming.
}
\author{
Guy J. Abel
}
\seealso{
\code{\link{ipf3.qi}}, \code{\link{fm}}
}
\examples{
## create P1 and P2 stock tables
dn <- LETTERS[1:4]
P1 <- matrix(c(1000, 100, 10, 0, 55, 555, 50, 5, 80, 40, 800, 40, 20, 25, 20, 200), 4, 4, 
        dimnames = list(pob = dn, por = dn), byrow = TRUE)
P2 <- matrix(c(950, 100, 60, 0, 80, 505, 75, 5, 90, 30, 800, 40, 40, 45, 0, 180), 4, 4, 
        dimnames = list(pob = dn, por = dn), byrow = TRUE)
# display with row and col totals
addmargins(P1)
addmargins(P2)

# no births and deaths
b <- rep(0, 4)
d <- rep(0, 4)

y <- ffs(P1, P2, d, b)
# display with row, col and table totals
round(addmargins(y$mu), 1)
# display with row and col totals
round(fm(y$mu), 1)

## alternative offset term
dis <- array(c(1, 2, 3, 4, 2, 1, 5, 6, 3, 4, 1, 7, 4, 6, 7, 1), c(4, 4, 4))
y <- ffs(P1, P2, d, b, dis)
# display with row, col and table totals
round(addmargins(y$mu), 1)
# display with row and col totals
round(fm(y$mu), 1)

## alternative P2 and changes in population from natural increase
P2 <- matrix(c(1060, 60, 10, 10, 45, 540, 40, 0, 70, 75, 770, 70, 30, 30, 20, 230), 4, 4, 
        dimnames = list(pob = dn, por = dn), byrow = TRUE)
# display with row and col totals
addmargins(P2)
b <- c(80, 20, 40, 60)
d <- c(70, 30, 50, 10)

y <- ffs(P1, P2, d, b)
# display with row, col and table totals
round(addmargins(y$mu), 1)
# display with row and col totals
round(fm(y$mu), 1)
}
