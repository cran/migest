#' Annual origin destination migration flows between Korean regions alongside selected geographic, economic and demographic variables.
#'
#' Origin-destination migration flows between 2012 and 2020 based on first level administrative regions.
#'
#' @format Data frame with 2,601 rows and 20 columns:
#' \describe{
#'   \item{orig}{Origin region}
#'   \item{dest}{Destination region}
#'   \item{year}{Year of flow}
#'   \item{flow}{Migration flow. Data obtained from KOSIS}
#'   \item{dist_cent}{Distance (in km) between geographic centroids, calculated from \code{geosphere::distm()}}
#'   \item{dist_min}{Minimum distance (in km) between regions, calculated from \code{sf::st_distance()}}
#'   \item{dist_pw}{Distance (in km) between population weighted centroids, calculated from \code{geosphere::distm()} using WorldPop estimates of 2020 regional population centroids}
#'   \item{contig}{Indicate if regions share a border}
#'   \item{orig_pop}{Population (in millions) of origin region. Data obtained from KOSIS.}
#'   \item{dest_pop}{Population (in millions) of destination region. Data obtained from KOSIS.}
#'   \item{orig_area}{Geographic area (in km^2) of origin region, calculated from \code{sf::st_area()}}
#'   \item{dest_area}{Geographic area (in km^2) of destination region, calculated from \code{sf::st_area()}}
#'   \item{orig_gdp_pc}{GDP per capita of origin region. Data obtained from KOSIS.}
#'   \item{orig_ginc_pc}{Gross regional income per capita of origin region. Data obtained from KOSIS.}
#'   \item{orig_iinc_pc}{Individual income per capita of origin region. Data obtained from KOSIS.}
#'   \item{orig_pconsum_pc}{Personal consumption per capita of origin region. Data obtained from KOSIS.}
#'   \item{dest_gdp_pc}{GDP per capita of destination region. Data obtained from KOSIS.}
#'   \item{dest_ginc_pc}{Gross regional income per capita of destination region. Data obtained from KOSIS.}
#'   \item{dest_iinc_pc}{Individual income per capita of destination region. Data obtained from KOSIS.}
#'   \item{dest_pconsum_pc}{Personal consumption per capita of destination region. Data obtained from KOSIS.}
#' }
#' @source Statistics Korea, Internal Migration Statistics. Data downloaded from \url{https://kosis.kr/eng} in July 2021.
#' @source Robin Edwards, Maksym Bondarenko, Andrew J. Tatem and Alessandro Sorichetta. Unconstrained subnational Population Weighted Density in 2000, 2005, 2010, 2015 and 2020 ( 100m resolution ). WorldPop, University of Southampton, UK. 
#' @source Source: Statistics Korea, Population Statistics Based on Resident Registration. Data downloaded from \url{https://kosis.kr/eng} in July 2021.
#' @source Source: Statistics Korea, Regional GDP, Gross regional income and Individual income. Data downloaded from \url{https://kosis.kr/eng} in November 2023.
"korea_gravity"