#' Prepare data for disaggregation modelling
#' 
#' \emph{prepare_data} function is used to extract all the data required for fitting a disaggregation model. 
#' Designed to be used in the \emph{disaggregation::fit_model} function.
#' 
#' Takes a SpatialPolygonDataFrame with the response data and a RasterStack of covariates. 
#' 
#' Extract the values of the covariates (as well as the aggregation raster, if given) at each pixel within the polygons 
#' (\emph{parallelExtract} function). This is done in parallel and \emph{n.cores} argument is used to set the number of cores 
#' to use for covariate extraction. This can be the number of covariates used in the model.
#' 
#' The aggregation raster defines how the pixels within each polygon are aggregated. 
#' The disaggregation model performs a weighted sum of the pixel prediction, weighted by the pixel values in the aggregation raster.
#' For disease incidence rate you use the population raster to aggregate pixel incidence rate by summing the number of cases 
#' (rate weighted by population). If no aggregation raster is provided a uniform distribution is assumed, i.e. the pixel predictions 
#' are aggregated to polygon level by summing the pixel values.
#' 
#' Makes a matrix that contains the start and end pixel index for each polygon. Builds an INLA mesh to use for the spatial field 
#' (\emph{getStartendindex} function).
#' 
#' The \emph{mesh.args} argument allows you to supply a list of INLA mesh parameters to control the mesh used for the spatial field 
#' (\emph{build_mesh} function).
#' 
#' The \emph{na.action} flag is automatically off. If there are any NAs in the response or covariate data within the polygons the 
#' \emph{prepare_data} method will error. Ideally the NAs in the data would be dealt with beforehand, however, setting na.action = TRUE 
#' will automatically deal with NAs. It removes any polygons that have NAs as a response, sets any aggregation pixels with NA to zero 
#' and sets covariate NAs pixels to the median value for the that covariate.
#' 
#' @param x sf object containing at least three columns: one with the geometried, one with the id for the polygons (\emph{id_var}) and one with the response count data (\emph{response_var}); for binomial data, i.e survey data, it can also contain a sample size column (\emph{sample_size_var}).
#' @param covariate_rasters RasterStack of covariate rasters to be used in the model.
#' @param aggregation_raster Raster to aggregate pixel level predictions to polygon level e.g. population to aggregate prevalence. If this is not supplied a uniform raster will be used.
#' @param id_var Name of column in SpatialPolygonDataFrame object with the polygon id.
#' @param response_var Name of column in SpatialPolygonDataFrame object with the response data.
#' @param sample_size_var For survey data, name of column in SpatialPolygonDataFrame object (if it exists) with the sample size data.
#' @param mesh.args list of parameters that control the mesh structure with the same names as used by INLA.
#' @param na.action logical. If TRUE, NAs in response will be removed, covariate NAs will be given the median value, aggregation NAs will be set to zero. Default FALSE (NAs in response or covariate data within the polygons will give errors).
#' @param makeMesh logical. If TRUE, build INLA mesh, takes some time. Default TRUE.
#' @param ncores Number of cores used to perform covariate extraction.
#'
#' @return A list is returned of class \code{disag_data}. 
#' The functions \emph{summary}, \emph{print} and \emph{plot} can be used on \code{disag_data}. 
#' The list  of class \code{disag_data} contains:
#'  \item{x }{The sf object used as an input.} 
#'  \item{covariate_rasters }{The SpatRaster used as an input.} 
#'  \item{polygon_data }{A data frame with columns of \emph{area_id}, \emph{response} and \emph{N} (sample size: all NAs unless using binomial data). Each row represents a polygon.}
#'  \item{covariate_data }{A data frame with columns of \emph{area_id}, \emph{cell_id} and one for each covariate in \emph{covariate_rasters}. Each row represents a pixel in a polygon.}
#'  \item{aggregation_pixels }{An array with the value of the aggregation raster for each pixel in the same order as the rows of \emph{covariate_data}.}
#'  \item{coordsForFit }{A matrix with two columns of x, y coordinates of pixels within the polygons. Used to make the spatial field.}
#'  \item{coordsForPrediction }{A matrix with two columns of x, y coordinates of pixels in the whole Raster. Used to make predictions.}
#'  \item{startendindex }{A matrix with two columns containing the start and end index of the pixels within each polygon.}
#'  \item{mesh }{A INLA mesh to be used for the spatial field of the disaggregation model.}
#' @import splancs
#' @import utils
#' @name prepare_data
#'
#' @examples 
#' \donttest{
#'  polygons <- list()
#'  for(i in 1:100) {
#'   row <- ceiling(i/10)
#'   col <- ifelse(i %% 10 != 0, i %% 10, 10)
#'   xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
#'   polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
#'  }
#' 
#'  polys <- do.call(raster::spPolygons, polygons)
#'  response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
#'  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
#' 
#'  r <- raster::raster(ncol=20, nrow=20)
#'  r <- raster::setExtent(r, raster::extent(spdf))
#'  r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
#'  r2 <- raster::raster(ncol=20, nrow=20)
#'  r2 <- raster::setExtent(r2, raster::extent(spdf))
#'  r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
#'  cov_rasters <- raster::stack(r, r2)
#' 
#'  test_data <- prepare_data(x = spdf, 
#'                            covariate_rasters = cov_rasters)
#' } 
#'                    
#' @export
#' 
#' 

prepare_data <- function(x, 
                         covariate_rasters,
                         aggregation_raster = NULL,
                         id_var = 'area_id', 
                         response_var = 'response', 
                         sample_size_var = NULL,
                         mesh.args = NULL, 
                         na.action = FALSE,
                         makeMesh = TRUE,
                         ncores = 2) {

  stopifnot(inherits(x, 'sf'))
  stopifnot(inherits(covariate_rasters, 'SpatRaster'))
  if(!is.null(aggregation_raster)) stopifnot(inherits(aggregation_raster, 'SpatRaster'))
  stopifnot(inherits(id_var, 'character'))
  stopifnot(inherits(response_var, 'character'))
  if(!is.null(mesh.args)) stopifnot(inherits(mesh.args, 'list'))
  
  # Check for NAs in response data
  na_rows <- is.na(x[, response_var, drop = TRUE])
  if(sum(na_rows) != 0) {
    if(na.action) {
      x <- x[!na_rows, ]
    } else {
      stop('There are NAs in the response data. Please deal with these, or set na.action = TRUE')
    }
  }
  
  polygon_data <- getPolygonData(x, id_var, response_var, sample_size_var)
  

  # Save raster layer names so we can reassign it to make sure names don't change.
  cov_names <- names(covariate_rasters)

  # If no aggregation raster is given, use a 'unity' raster
  if(is.null(aggregation_raster)) {
    aggregation_raster <- covariate_rasters[[1]]
    terra::values(aggregation_raster) <- rep(1, terra::ncell(aggregation_raster))
  }
  names(aggregation_raster) <- 'aggregation_raster'

  
  covariate_rasters <- c(covariate_rasters, aggregation_raster)
  covariate_data <- parallelExtract(covariate_rasters, x, fun = NULL, id = id_var)

  # Remove the aggregation raster
  covariate_rasters <- covariate_rasters[[seq(nlyr(covariate_rasters) - 1)]]
  
  names(covariate_rasters) <- cov_names
  
  aggregation_pixels <- as.numeric(covariate_data[ , ncol(covariate_data)])
  covariate_data <- covariate_data[, -ncol(covariate_data)]
  
  # Check for NAs in population data
  if(sum(is.na(aggregation_pixels)) != 0) {
    if(na.action) {
      aggregation_pixels[is.na(aggregation_pixels)] <- 0
    } else {
      stop('There are NAs in the aggregation rasters within polygons. Please deal with these, or set na.action = TRUE')
    }
  }
  
  # Check for NAs in covariate data
  if(sum(is.na(covariate_data)) != 0) {
    if(na.action) {
      covariate_data[-c(1:2)] <- sapply(covariate_data[-c(1:2)], function(x) { x[is.na(x)] <- stats::median(x, na.rm = T); return(x) })
    } else {
      stop('There are NAs in the covariate rasters within polygons. Please deal with these, or set na.action = TRUE')
    }
  }
  
  coordsForFit <- extractCoordsForMesh(covariate_rasters, selectIds = covariate_data$cellid)
  
  coordsForPrediction <- extractCoordsForMesh(covariate_rasters)
  
  startendindex <- getStartendindex(covariate_data, polygon_data, id_var = id_var)
  
  if(makeMesh) {
    if(!requireNamespace('INLA', quietly = TRUE)) {
      mesh <- NULL
      message("Cannot build mesh as INLA is not installed. If you need a spatial field in your model, you must install INLA.")
    } else {
      mesh <- build_mesh(x, mesh.args)
    }
  } else {
    mesh <- NULL
    message("A mesh is not being built. You will not be able to run a spatial model without a mesh.")
  }
  
  disag_data <- list(x = x,
                     shapefile_names = list(id_var = id_var, response_var = response_var),
                     covariate_rasters = covariate_rasters,
                     polygon_data = polygon_data,
                     covariate_data = covariate_data,
                     aggregation_pixels = aggregation_pixels,
                     coordsForFit = coordsForFit,
                     coordsForPrediction = coordsForPrediction,
                     startendindex = startendindex,
                     mesh = mesh)
  
  class(disag_data) <- c('disag_data', 'list')
  
  return(disag_data)
  
}

#' Function to fit the disaggregation model
#'
#' @param x SpatialPolygonDataFrame containing the response data 
#' @param shapefile_names List of 2: polygon id variable name and response variable name from x
#' @param covariate_rasters RasterStack of covariates
#' @param polygon_data data.frame with two columns: polygon id and response
#' @param covariate_data data.frame with cell id, polygon id and covariate columns
#' @param aggregation_pixels vector with value of aggregation raster at each pixel
#' @param coordsForFit coordinates of the covariate data points within the polygons in x
#' @param coordsForPrediction coordinates of the covariate data points in the whole raster extent
#' @param startendindex matrix containing the start and end index for each polygon
#' @param mesh inla.mesh object to use in the fit
#' 
#' @return A list is returned of class \code{disag_data}. 
#' The functions \emph{summary}, \emph{print} and \emph{plot} can be used on \code{disag_data}. 
#' The list  of class \code{disag_data} contains:
#'  \item{x }{The SpatialPolygonDataFrame used as an input.} 
#'  \item{covariate_rasters }{The RasterStack used as an input.} 
#'  \item{polygon_data }{A data frame with columns of \emph{area_id}, \emph{response} and \emph{N} (sample size: all NAs unless using binomial data). Each row represents a polygon.}
#'  \item{covariate_data }{A data frame with columns of \emph{area_id}, \emph{cell_id} and one for each covariate in \emph{covariate_rasters}. Each row represents a pixel in a polygon.}
#'  \item{aggregation_pixels }{An array with the value of the aggregation raster for each pixel in the same order as the rows of \emph{covariate_data}.}
#'  \item{coordsForFit }{A matrix with two columns of x, y coordinates of pixels within the polygons. Used to make the spatial field.}
#'  \item{coordsForPrediction }{A matrix with two columns of x, y coordinates of pixels in the whole Raster. Used to make predictions.}
#'  \item{startendindex }{A matrix with two columns containing the start and end index of the pixels within each polygon.}
#'  \item{mesh }{A INLA mesh to be used for the spatial field of the disaggregation model.}
#'
#' @name as.disag_data
#' 
#' @export


as.disag_data <- function(x, 
                          shapefile_names,
                          covariate_rasters, 
                          polygon_data, 
                          covariate_data, 
                          aggregation_pixels,
                          coordsForFit, 
                          coordsForPrediction,
                          startendindex, 
                          mesh = NULL) {
  
  stopifnot(inherits(x, 'SpatialPolygonsDataFrame'))
  stopifnot(inherits(shapefile_names, 'list'))
  stopifnot(inherits(covariate_rasters, c('RasterBrick', 'RasterStack')))
  stopifnot(inherits(polygon_data, 'data.frame'))
  stopifnot(inherits(covariate_data, 'data.frame'))
  stopifnot(inherits(aggregation_pixels, 'numeric'))
  stopifnot(inherits(coordsForFit, 'matrix'))
  stopifnot(inherits(coordsForPrediction, 'matrix'))
  stopifnot(inherits(startendindex, 'matrix'))
  if(!is.null(mesh)) {
    stopifnot(inherits(mesh, 'inla.mesh'))
  }
  
  disag_data <- list(x = x,
                     shapefile_names = shapefile_names,
                     covariate_rasters = covariate_rasters,
                     polygon_data = polygon_data,
                     covariate_data = covariate_data,
                     aggregation_pixels = aggregation_pixels,
                     coordsForFit = coordsForFit,
                     coordsForPrediction = coordsForPrediction,
                     startendindex = startendindex,
                     mesh = mesh)
  
  class(disag_data) <- c('disag_data', 'list')
  
  return(disag_data)
}
