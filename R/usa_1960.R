#' US population totals in 1950 and 1960 by place of birth, age, sex and race
#'
#' Population data by place of birth, age, sex and race in 1950 and 1960 
#'
#' @format Data frame with 288 rows and 7 columns:
#' \describe{
#'   \item{birthplace}{Place of birth (US Census area)}
#'   \item{race}{Race from `white` or `non-white`}
#'   \item{sex}{Sex from `male` or `female`}
#'   \item{age_1950}{Age group in 1950}
#'   \item{age_1960}{Age group in 1960}
#'   \item{pop_1950}{Enumerated population in 1950}
#'   \item{pop_1960}{Enumerated population in 1960}
#' }
#' @source Data scraped from Table D, pp. 183-191 of Eldridge, H., & Kim, Y. (1968). The estimation of intercensal migration from birth-residence statistics: a study of data for the United States, 1950 and 1960 (PSC Analytical and Technical Report Series, Issue 7). \url{https://repository.upenn.edu/entities/publication/2a11a5f7-3ddf-47f3-a47d-1de5254f4cc5}
"usa_1960"