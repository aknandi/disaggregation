#' Prepare data for disaggregation modelling
#'
#' \emph{prepare_data} function is used to extract all the data required for fitting a disaggregation model.
#' Designed to be used in the \emph{disaggregation::disag_model} function.
#'
#' Takes a sf object with the response data and a SpatRaster of covariates.
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
#' @param polygon_shapefile sf object containing at least three columns: one with the geometried, one with the id for the polygons (\emph{id_var}) and one with the response count data (\emph{response_var}); for binomial data, i.e survey data, it can also contain a sample size column (\emph{sample_size_var}).
#' @param covariate_rasters SpatRaster of covariate rasters to be used in the model.
#' @param aggregation_raster SpatRaster to aggregate pixel level predictions to polygon level e.g. population to aggregate prevalence. If this is not supplied a uniform raster will be used.
#' @param id_var Name of column in sf object with the polygon id.
#' @param response_var Name of column in sf object with the response data.
#' @param sample_size_var For survey data, name of column in sf object (if it exists) with the sample size data.
#' @param mesh_args list of parameters that control the mesh structure with the same names as used by INLA.
#' @param na_action logical. If TRUE, NAs in response will be removed, covariate NAs will be given the median value, aggregation NAs will be set to zero. Default FALSE (NAs in response or covariate data within the polygons will give errors).
#' @param make_mesh logical. If TRUE, build INLA mesh, takes some time. Default TRUE.
#' @param mesh.args Deprecated.
#' @param na.action Deprecated.
#' @param makeMesh Deprecated.
#' @param ncores Deprecated.
#'
#' @return A list is returned of class \code{disag_data}.
#' The functions \emph{summary}, \emph{print} and \emph{plot} can be used on \code{disag_data}.
#' The list  of class \code{disag_data} contains:
#'  \item{x }{The sf object used as an input.}
#'  \item{covariate_rasters }{The SpatRaster used as an input.}
#'  \item{polygon_data }{A data frame with columns of \emph{area_id}, \emph{response} and \emph{N} (sample size: all NAs unless using binomial data). Each row represents a polygon.}
#'  \item{covariate_data }{A data frame with columns of \emph{area_id}, \emph{cell_id} and one for each covariate in \emph{covariate_rasters}. Each row represents a pixel in a polygon.}
#'  \item{aggregation_pixels }{An array with the value of the aggregation raster for each pixel in the same order as the rows of \emph{covariate_data}.}
#'  \item{coords_for_fit }{A matrix with two columns of x, y coordinates of pixels within the polygons. Used to make the spatial field.}
#'  \item{coords_for_prediction }{A matrix with two columns of x, y coordinates of pixels in the whole Raster. Used to make predictions.}
#'  \item{start_end_index }{A matrix with two columns containing the start and end index of the pixels within each polygon.}
#'  \item{mesh }{A INLA mesh to be used for the spatial field of the disaggregation model.}
#' @import splancs
#' @import utils
#' @name prepare_data
#'
#' @examples
#' \donttest{
#' polygons <- list()
#' for(i in 1:100) {
#'   row <- ceiling(i/10)
#'   col <- ifelse(i %% 10 != 0, i %% 10, 10)
#'   xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
#'   polygons[[i]] <- list(cbind(c(xmin, xmax, xmax, xmin, xmin),
#'                               c(ymax, ymax, ymin, ymin, ymax)))
#' }
#'
#' polys <- lapply(polygons,sf::st_polygon)
#' response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
#' spdf <- sf::st_sf(response_df,geometry=polys)
#'
#' r <- terra::rast(nrow=20,ncol=20)
#' terra::ext(r) <- terra::ext(spdf)
#' r[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
#'
#' r2 <- terra::rast(nrow=20,ncol=20)
#' terra::ext(r2) <- terra::ext(spdf)
#' r2[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
#' cov_rasters <- c(r, r2)
#'
#' test_data <- prepare_data(polygon_shapefile = spdf,
#'                           covariate_rasters = cov_rasters)
#' }
#'
#' @export
#'
#'

prepare_data <- function(polygon_shapefile,
                         covariate_rasters,
                         aggregation_raster = NULL,
                         id_var = 'area_id',
                         response_var = 'response',
                         sample_size_var = NULL,
                         mesh_args = NULL,
                         na_action = FALSE,
                         make_mesh = TRUE,
                         mesh.args = NULL,
                         na.action = NULL,
                         makeMesh = NULL,
                         ncores = NULL) {

  # Deal with deprecated parameters

  if (!is.null(na.action) && missing(na_action)) {
    na_action <- na.action
    message("na.action is deprecated and will be removed in a future version - please use na_action instead")
  }

  if (!is.null(mesh.args) && missing(mesh_args)) {
    mesh_args <- mesh.args
    message("mesh.args is deprecated and will be removed in a future version - please use mesh_args instead")
  }

  if (!is.null(makeMesh) && missing(make_mesh)) {
    make_mesh <- makeMesh
    message("makeMesh is deprecated and will be removed in a future version - please use make_mesh instead")
  }

  if (!missing("ncores"))
    warning("ncores is deprecated and will be removed in a future version")

  stopifnot(inherits(polygon_shapefile, 'sf'))
  stopifnot(inherits(covariate_rasters, 'SpatRaster'))
  if(!is.null(aggregation_raster)) stopifnot(inherits(aggregation_raster, 'SpatRaster'))
  stopifnot(inherits(id_var, 'character'))
  stopifnot(inherits(response_var, 'character'))
  if(!is.null(mesh_args)) stopifnot(inherits(mesh_args, 'list'))

  # Check for NAs in response data
  na_rows <- is.na(polygon_shapefile[, response_var, drop = TRUE])
  if(sum(na_rows) != 0) {
    if(na_action) {
      polygon_shapefile <- polygon_shapefile[!na_rows, ]
    } else {
      stop('There are NAs in the response data. Please deal with these, or set na.action = TRUE')
    }
  }

  # If no aggregation raster is given, use a 'unity' raster
  if(is.null(aggregation_raster)) {
    aggregation_raster <- covariate_rasters[[1]]
    terra::values(aggregation_raster) <- rep(1, terra::ncell(aggregation_raster))
  }
  names(aggregation_raster) <- 'aggregation_raster'

  # Check for polygons with zero aggregation
  aggregation_sum <- terra::extract(aggregation_raster, polygon_shapefile, cells = TRUE, na.rm = TRUE, ID = TRUE, fun = "sum")
  zero_aggregation <- aggregation_sum[aggregation_sum[, 2] == 0,]
  if (nrow(zero_aggregation) > 0){
    zero_polygons <- polygon_shapefile[zero_aggregation$ID,]
    if (na_action){
      for (p in zero_polygons[[id_var]]){
        message(paste0(p, " has zero aggregation values and has been removed from the response data"))
      }
      polygon_shapefile <- polygon_shapefile[-zero_aggregation$ID,]
    } else {
      for (p in zero_polygons[[id_var]]){
        message(paste0(p, " has zero aggregation values"))
      }
      stop("Please remove the polygons from the response data or set na.action = TRUE")
    }
  }

  polygon_data <- getPolygonData(polygon_shapefile, id_var, response_var, sample_size_var)

  # Check for non-numeric covariate values
  cv_df <- terra::as.data.frame(covariate_rasters, xy = FALSE)
  cv_classes <- unlist(lapply(cv_df, class))
  cv_check <- all(as.vector(cv_classes) == "numeric")
  if (!cv_check){
    non_numeric <- which(cv_classes != "numeric")
    for (raster in non_numeric){
      warning(paste0("The values of ", names(covariate_rasters)[raster], " are not numeric"))
    }
  }

  # Save raster layer names so we can reassign it to make sure names don't change.
  cov_names <- names(covariate_rasters)

  covariate_rasters <- c(covariate_rasters, aggregation_raster)
  covariate_data <- terra::extract(covariate_rasters, terra::vect(polygon_shapefile), cells=TRUE, na.rm=TRUE, ID=TRUE)
  #merge to transfer area_id and then tidy up
  polygon_data$area_n <- 1:nrow(polygon_data)
  covariate_data <- merge(covariate_data, polygon_data, by.x = "ID", by.y = "area_n")
  covariate_data <- covariate_data[ , !(names(covariate_data) %in% c("ID", "response", "N"))]
  colnames(covariate_data )[colnames(covariate_data ) == "area_id"] <- id_var
  polygon_data <- polygon_data[ , !(names(polygon_data) %in% c("area_n"))]

  # Remove the aggregation raster
  cov_filter <- !(names(covariate_data) %in% c('aggregation_raster'))
  covariate_rasters <- covariate_rasters[[cov_filter]]
  names(covariate_rasters) <- cov_names

  agg_filter <- names(covariate_data) %in% c('aggregation_raster')
  aggregation_pixels <- as.numeric(covariate_data[ , agg_filter])
  covariate_data <- covariate_data[, !agg_filter]

  # Check for NAs in population data
  if(sum(is.na(aggregation_pixels)) != 0) {
    if(na_action) {
      aggregation_pixels[is.na(aggregation_pixels)] <- 0
    } else {
      stop('There are NAs in the aggregation rasters within polygons. Please deal with these, or set na_action = TRUE')
    }
  }

  # Check for NAs in covariate data
  if(sum(is.na(covariate_data)) != 0) {
    if(na_action) {
      cov_filter <- !(names(covariate_data) %in% c(id_var,'cell'))
      covariate_data[ , cov_filter] <- sapply(covariate_data[ , cov_filter], function(x) { x[is.na(x)] <- stats::median(x, na.rm = T); return(x) })
    } else {
      stop('There are NAs in the covariate rasters within polygons. Please deal with these, or set na_action = TRUE')
    }
  }

  coords_for_fit <- extractCoordsForMesh(covariate_rasters, selectIds = covariate_data$cell)

  coords_for_prediction <- extractCoordsForMesh(covariate_rasters)

  start_end_index <- getStartendindex(covariate_data, polygon_data, id_var = id_var)

  if(make_mesh) {
      mesh <- build_mesh(polygon_shapefile, mesh_args)
    } else {
    mesh <- NULL
    message("A mesh is not being built. You will not be able to run a spatial model without a mesh.")
  }

  disag_data <- list(polygon_shapefile = polygon_shapefile,
                     shapefile_names = list(id_var = id_var, response_var = response_var),
                     covariate_rasters = covariate_rasters,
                     polygon_data = polygon_data,
                     covariate_data = covariate_data,
                     aggregation_pixels = aggregation_pixels,
                     coords_for_fit = coords_for_fit,
                     coords_for_prediction = coords_for_prediction,
                     start_end_index = start_end_index,
                     mesh = mesh)

  class(disag_data) <- c('disag_data', 'list')

  return(disag_data)

}

#' Function to fit the disaggregation model
#'
#' @param polygon_shapefile sf object containing the response data
#' @param shapefile_names List of 2: polygon id variable name and response variable name from x
#' @param covariate_rasters SpatRaster of covariates
#' @param polygon_data data.frame with two columns: polygon id and response
#' @param covariate_data data.frame with cell id, polygon id and covariate columns
#' @param aggregation_pixels vector with value of aggregation raster at each pixel
#' @param coords_for_fit coordinates of the covariate data points within the polygons in x
#' @param coords_for_prediction coordinates of the covariate data points in the whole raster extent
#' @param start_end_index matrix containing the start and end index for each polygon
#' @param coordsForFit Deprecated.
#' @param coordsForPrediction Deprecated.
#' @param startendindex Deprecated.
#' @param mesh inla.mesh object to use in the fit
#'
#' @return A list is returned of class \code{disag_data}.
#' The functions \emph{summary}, \emph{print} and \emph{plot} can be used on \code{disag_data}.
#' The list  of class \code{disag_data} contains:
#'  \item{x }{The sf object used as an input.}
#'  \item{covariate_rasters }{The SpatRaster used as an input.}
#'  \item{polygon_data }{A data frame with columns of \emph{area_id}, \emph{response} and \emph{N} (sample size: all NAs unless using binomial data). Each row represents a polygon.}
#'  \item{covariate_data }{A data frame with columns of \emph{area_id}, \emph{cell_id} and one for each covariate in \emph{covariate_rasters}. Each row represents a pixel in a polygon.}
#'  \item{aggregation_pixels }{An array with the value of the aggregation raster for each pixel in the same order as the rows of \emph{covariate_data}.}
#'  \item{coords_for_fit }{A matrix with two columns of x, y coordinates of pixels within the polygons. Used to make the spatial field.}
#'  \item{coords_for_prediction }{A matrix with two columns of x, y coordinates of pixels in the whole Raster. Used to make predictions.}
#'  \item{start_end_index }{A matrix with two columns containing the start and end index of the pixels within each polygon.}
#'  \item{mesh }{A INLA mesh to be used for the spatial field of the disaggregation model.}
#'
#' @name as.disag_data
#'
#' @export


as.disag_data <- function(polygon_shapefile,
                          shapefile_names,
                          covariate_rasters,
                          polygon_data,
                          covariate_data,
                          aggregation_pixels,
                          coords_for_fit,
                          coords_for_prediction,
                          start_end_index,
                          mesh = NULL,
                          coordsForFit = NULL,
                          coordsForPrediction = NULL,
                          startendindex = NULL) {

  # Handle deprecated variables
  if (!is.null(coordsForFit) && missing(coords_for_fit)) {
    coords_for_fit <- coordsForFit
    message("coordsForFit is deprecated and will be removed in a future version - please use coords_for_fit instead")
  }

  if (!is.null(coordsForPrediction) && missing(coords_for_prediction)) {
    coords_for_prediction <- coordsForPrediction
    message("coordsForPrediction is deprecated and will be removed in a future version - please use coords_for_prediction instead")
  }

  if (!is.null(startendindex) && missing(start_end_index)) {
    start_end_index <- startendindex
    message("startendindex is deprecated and will be removed in a future version - please use start_end_index instead")
  }

  stopifnot(inherits(polygon_shapefile, 'sf'))
  stopifnot(inherits(shapefile_names, 'list'))
  stopifnot(inherits(covariate_rasters, 'SpatRaster'))
  stopifnot(inherits(polygon_data, 'data.frame'))
  stopifnot(inherits(covariate_data, 'data.frame'))
  stopifnot(inherits(aggregation_pixels, 'numeric'))
  stopifnot(inherits(coords_for_fit, 'matrix'))
  stopifnot(inherits(coords_for_prediction, 'matrix'))
  stopifnot(inherits(start_end_index, 'matrix'))
  if(!is.null(mesh)) {
    stopifnot(inherits(mesh, 'inla.mesh'))
  }

  disag_data <- list(polygon_shapefile = polygon_shapefile,
                     shapefile_names = shapefile_names,
                     covariate_rasters = covariate_rasters,
                     polygon_data = polygon_data,
                     covariate_data = covariate_data,
                     aggregation_pixels = aggregation_pixels,
                     coords_for_fit = coords_for_fit,
                     coords_for_prediction = coords_for_prediction,
                     start_end_index = start_end_index,
                     mesh = mesh)

  class(disag_data) <- c('disag_data', 'list')

  return(disag_data)
}
