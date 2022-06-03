#' Lucas County House Price Data
#'
#' A dataset containing the prices and other attributes of 25,357 houses in
#' Lucas County, Ohio. The selling dates span years 1993 to 1998. Data taken
#' from \code{\link[spData]{house}} (\code{spData} package) and slightly modified to a \code{data.frame}.
#'
#' @format A data frame with 25357 rows and 25 variables:
#' \describe{
#'   \item{price}{(\code{integer}) selling price, in US dollars}
#'   \item{yrbuilt}{(\code{integer}) year the house was built}
#'   \item{stories}{(\code{factor}) levels are \code{"one", "bilevel",
#'   "multilvl",  "one+half", "two", "two+half", "three"}}
#'   \item{TLA}{(\code{integer}) total living area, in square feet.}
#'   \item{wall}{(\code{factor}) levels are \code{"stucdrvt", "ccbtile",
#'   "metlvnyl", "brick", "stone", "wood", "partbrk"}}
#'   \item{beds, baths, halfbaths}{(\code{integer}) number of corresponding
#'   rooms / facilities.}
#'   \item{frontage, depth}{dimensions of the lot. Unit is feet.}
#'   \item{garage}{(\code{factor}) levels are \code{"no garage", "basement",
#'   "attached", "detached", "carport"}}
#'   \item{garagesqft}{(\code{integer}) garage area, in square feet. If
#'   \code{garage == "no garage"}, then \code{garagesqft == 0}.}
#'   \item{rooms}{(\code{integer}) number of rooms}
#'   \item{lotsize}{(\code{integer}) area of lot, in square feet}
#'   \item{sdate}{(\code{Date}) selling date, in format \code{yyyy-mm-dd}}
#'   \item{avalue}{(\code{int}) appraised value}
#'   \item{s1993, s1994, s1995, s1996, s1997, s1998}{(\code{int}) dummies for
#'   selling year.}
#'   \item{syear}{(\code{factor}) levels are selling years \code{"1993", "1994",
#'   "1995", "1996", "1997", "1998"}}
#'   \item{long, lat}{(\code{numeric}) location of houses. Longitude and
#'   Latitude are given in \code{CRS(+init=epsg:2834)}, the Ohio North State
#'   Plane. Units are meters.}
#' }
#' @source \url{http://www.spatial-econometrics.com/html/jplv6.zip}
"house"


#' Sampled SVC Data
#'
#' A list object that contains sampled data of 500 observations. The data has 
#' been sampled using the \code{RandomFields} package. It is given in the list object 
#' \code{SVCdata} which contains the following.
#'
#' @format A `list` with the following entries:
#' \describe{
#'   \item{y}{(\code{numeric}) Response}
#'   \item{X}{(\code{numeric}) Covariates; first columns contains ones to model
#'   an intercept, the second column contains standard-normal sampled data.}
#'   \item{beta}{(\code{numeric}) The sampled Gaussian processes, which are
#'   usually unobserved. It uses a Matern covariance function and the true 
#'   parameters are given in the entry `true_pars`.}
#'   \item{eps}{(\code{numeric}) Error (or Nugget effect), i.e., drawn from a 
#'   zero-mean normal distribution with 0.5 standard deviation.}
#'   \item{locs}{(\code{numeric}) Locations sampled from a uniform distribution
#'   on the interval 0 to 10.}
#'   \item{true_pars}{(\code{data.frame}) True parameters of the GP-based SVC 
#'   model with Gaussian process mean, variance, and range. Additionally, the 
#'   smoothness (nu) is given.}
#' }
"SVCdata"
