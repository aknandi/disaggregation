#' Function to match pixels to their corresponding polygon
#' 
#' From the covariate data and polygon data, the function matches the polygon id between the two to find 
#' which pixels from the covariate data are contained in each of the polygons.
#' 
#' Takes a data.frame containing the covariate data with a polygon id column and one column for each covariate, 
#' and another data.frame containing polygon data with a polygon id, response and sample size column (as returned 
#' by \code{getPolygonData} function).
#' 
#' Returns a matrix with two columns and one row for each polygon. The first column is the index of the first row in
#' covariate data that corresponds to that polygon, the second column is the index of the last row in
#' covariate data that corresponds to that polygon.
#'
#' @param covariates data.frame with each covariate as a column an and id column.
#' @param polygon_data data.frame with polygon id and response data.
#' @param id_var string with the name of the column in the covariate data.frame containing the polygon id.
#' 
#' @return A matrix with two columns and one row for each polygon. The first column is the index of the first row in
#' covariate data that corresponds to that polygon, the second column is the index of the last row in
#' covariate data that corresponds to that polygon.
#'
#' @name getStartendindex
#'
#' @examples {
#'  covs <- data.frame(area_id = c(1, 1, 1, 2, 2, 3, 3, 3, 3), response = c(3, 9, 5, 2, 3, 6, 7, 3, 5))
#'  response <- data.frame(area_id = c(1, 2, 3), response = c(4, 7, 2), N = c(NA, NA, NA))
#'  getStartendindex(covs, response, 'area_id')
#' }
#'
#'
#' @export

getStartendindex <- function(covariates, polygon_data, id_var = 'area_id') {

  stopifnot(ncol(polygon_data) == 3)
  stopifnot(ncol(covariates) >= 2)
  stopifnot(nrow(covariates) > nrow(polygon_data))
  stopifnot(sum(polygon_data$area_id %in% covariates[, id_var]) == nrow(polygon_data))

  # Create  startendindex matrix
  # This defines which pixels in the matrix are associated with which polygon.
  startendindex <- lapply(unique(covariates[, id_var]), function(x) range(which(covariates[, id_var] == x)))

  startendindex <- do.call(rbind, startendindex)

  whichindices <- match(polygon_data$area_id, unique(covariates[, id_var]))

  # c++ is zero indexed.
  startendindex <- startendindex[whichindices, ] - 1L

  return(startendindex)
}

