#' @title Extract Number of Unique Locations
#'
#' @description Function to extract the number of unique locations in the data 
#' set used in an MLE of the \code{\link{SVC_mle}} object.
#'
#' @param object \code{\link{SVC_mle}} object
#'
#' @return integer with the number of unique locations
#'
#' @author Jakob Dambon
#'
#'
#' @export
nlocs <- function(object) {
  nrow(unique(object$MLE$call.args$locs))
}

