#' Find the range of pixels in each polygon
#'
#' @param covariates data.frame with each covariate as a column an and id column
#' @param polygon_data data.frame with polygon id and response data
#' @param id_var string with the name of the column containing the polygon id
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

