#' Extract polygon id and response data into a data.frame from a SpatialPolygonsDataFrame
#'
#' Returns a data.frame with a row for each polygon in the SpatialPolygonDataFrame and columns: area_id, response and N, containing the id of the
#' polygon, the values of the response for that polygon, and the sample size respectively. If the data is not survey data (the sample size does
#' not exist), this column will contain NAs.
#'
#' @param shape A sf object containing response data.
#' @param id_var Name of column in shape object with the polygon id. Default 'area_id'.
#' @param response_var Name of column in shape object with the response data. Default 'response'.
#' @param sample_size_var For survey data, name of column in SpatialPolygonDataFrame object (if it exists) with the sample size data. Default NULL.
#'
#' @return A data.frame with a row for each polygon in the SpatialPolygonDataFrame and columns: area_id, response and N, containing the id of the
#' polygon, the values of the response for that polygon, and the sample size respectively. If the data is not survey data (the sample size does
#' not exist), this column will contain NAs.
#'
#' @export
#' @examples {
#'  polygons <- list()
#'  for(i in 1:100) {
#'    row <- ceiling(i/10)
#'    col <- ifelse(i %% 10 != 0, i %% 10, 10)
#'    xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
#'    polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
#'  }
#'
#'  polys <- do.call(raster::spPolygons, polygons)
#'  response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
#'  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
#'
#'  getPolygonData(spdf, id_var = 'area_id', response_var = 'response')
#' }
#'
#'

getPolygonData <- function(shape, id_var = 'area_id', response_var = 'response', sample_size_var = NULL) {

  if(is.null(sample_size_var)) {
    polygon_df <- shape[, c(id_var, response_var), drop = TRUE]
    polygon_df$N <- rep(NA, nrow(polygon_df))
  } else {
    polygon_df <- shape[, c(id_var, response_var, sample_size_var), drop = TRUE]
  }

  names(polygon_df) <- c('area_id', 'response', 'N')

  return(polygon_df)
}


#' Get a RasterStack of covariates from a folder containing .tif files
#'
#' Looks in a specified folder for raster files. Returns a RasterStack of the rasters cropped to the extent specified by the shape parameter.
#'
#' @param directory Filepath to the directory containing the rasters.
#' @param file_pattern Pattern the filenames must match. Default is all files ending in .tif .
#' @param shape An object with an extent that the rasters will be cropped to.
#'
#' @return A RasterStack of the raster files in the directory
#'
#' @export
#' @examples
#' \dontrun{
#'   getCovariateRasters('/home/rasters', '.tif$', shape)
#'  }
#'

getCovariateRasters <- function(directory, file_pattern = '.tif$', shape) {

  stopifnot(dir.exists(directory))

  covariate_files <- list.files(directory, pattern = file_pattern, full.names = TRUE)
  stopifnot(length(covariate_files) != 0)

  covariate_rasters <- lapply(covariate_files, function(x) raster::raster(x))
  covariate_stack <- raster::stack(covariate_rasters)

  covariate_stack <- raster::crop(covariate_stack, shape)

  return(covariate_stack)
}

# Extract coordinates from raster to use constructing the INLA mesh
#
# @param cov_rasters RasterStack of the covariate rasters.
# @param selectIds numeric vector containing cell ids to retain. Default NULL retains all cell ids in the covariate rasters.
#
# @return A matrix containing the coordinates used to make the mesh

extractCoordsForMesh <- function(cov_rasters, selectIds = NULL) {

  stopifnot(inherits(cov_rasters, 'SpatRaster'))
  if(!is.null(selectIds)) stopifnot(inherits(selectIds, 'numeric'))

  points_raster <- cov_rasters[[1]]
  points_raster[is.na(terra::values(points_raster, mat = FALSE))] <- -9999
  raster_pts <- terra::as.points(points_raster)
  coords <- terra::crds(raster_pts)

  # If specified, only retain certain pixel ids
  if(!is.null(selectIds)) {
    coords <- coords[selectIds, ]
  }

  return(coords)

}



