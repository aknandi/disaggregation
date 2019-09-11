#' Function to collect all data together for disaggregation model
#'
#' @param polygon_shapefile SpatialPolygonDataFrame containing the response data 
#' @param covariate_rasters RasterStack of covariates
#' @param aggregation_raster Raster to aggregate predicted pixel values e.g. population to aggregate prevalence
#' @param id_var Name of column in SpatialPolygonDataFrame object with the polygon id
#' @param response_var Name of column in SpatialPolygonDataFrame object with the response data
#' @param mesh.args list of parameters that control the mesh structure with the same names as used by INLA
#' @param ncores Number of cores used to perform covariate extraction
#'
#' @name prepare_data
#'
#' @examples 
#' \dontrun{
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
#'  r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
#'  r2 <- raster::raster(ncol=20, nrow=20)
#'  r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
#'  cov_rasters <- raster::stack(r, r2)
#' 
#'  test_data <- prepare_data(polygon_shapefile = spdf, 
#'                            covariate_rasters = cov_rasters)
#' } 
#'                    
#' @export

prepare_data <- function(polygon_shapefile, 
                         covariate_rasters,
                         aggregation_raster = NULL,
                         id_var = 'area_id', 
                         response_var = 'response', 
                         mesh.args = NULL, 
                         ncores = 2) {

  stopifnot(inherits(polygon_shapefile, 'SpatialPolygonsDataFrame'))
  stopifnot(inherits(covariate_rasters, 'Raster'))
  if(!is.null(aggregation_raster)) stopifnot(inherits(aggregation_raster, 'Raster'))
  stopifnot(inherits(id_var, 'character'))
  stopifnot(inherits(response_var, 'character'))
  if(!is.null(mesh.args)) stopifnot(inherits(mesh.args, 'list'))
  
  polygon_data <- getPolygonData(polygon_shapefile, id_var = id_var, response_var = response_var)
  
  # If no aggregation raster is given, use a 'unity' raster
  if(is.null(aggregation_raster)) {
    aggregation_raster <- covariate_rasters[[1]]
    aggregation_raster <- raster::setValues(aggregation_raster, rep(1, raster::ncell(aggregation_raster)))
    names(aggregation_raster) <- 'aggregation_raster'
  }
  
  covariate_rasters <- raster::addLayer(covariate_rasters, aggregation_raster)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  covariate_data <- parallelExtract(covariate_rasters, polygon_shapefile, fun = NULL, id = id_var)
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  covariate_rasters <- raster::dropLayer(covariate_rasters, raster::nlayers(covariate_rasters))
  
  aggregation_pixels <- as.numeric(covariate_data[ , ncol(covariate_data)])
  covariate_data <- covariate_data[, -ncol(covariate_data)]
  
  coords <- extractCoordsForMesh(covariate_rasters, covariate_data)
  
  startendindex <- getStartendindex(covariate_data, polygon_data, id_var = id_var)
  
  mesh <- build_mesh(polygon_shapefile, mesh.args)
  
  disag_data <- list(polygon_shapefile = polygon_shapefile,
                     covariate_rasters = covariate_rasters,
                     polygon_data = polygon_data,
                     covariate_data = covariate_data,
                     aggregation_pixels = aggregation_pixels,
                     coords = coords,
                     startendindex = startendindex,
                     mesh = mesh)
  
  class(disag_data) <- c('disag.data', 'list')
  
  return(disag_data)
  
}

#' Function to fit the disaggregation model
#'
#' @param polygon_shapefile SpatialPolygonDataFrame containing the response data 
#' @param covariate_rasters RasterStack of covariates
#' @param polygon_data data.frame with two columns: polygon id and response
#' @param covariate_data data.frame with cell id, polygon id and covariate columns
#' @param aggregation_pixels vector with value of aggregation raster at each pixel
#' @param coords coordinates of the covariate data points
#' @param startendindex matrix containing the start and end index for each polygon
#' @param mesh inla.mesh object to use in the fit
#' 
#' @name as.disag.data
#'
#' @examples 
#' \dontrun{
#'   as.disag.data(polygon_data, covariate_data, coords, startendindex, mesh)
#'  }
#' 
#' @export


as.disag.data <- function(polygon_shapefile, 
                          covariate_rasters, 
                          polygon_data, 
                          covariate_data, 
                          aggregation_pixels,
                          coords, 
                          startendindex, 
                          mesh) {
  
  stopifnot(inherits(polygon_shapefile, 'SpatialPolygonsDataFrame'))
  stopifnot(inherits(covariate_rasters, c('RasterBrick', 'RasterStack')))
  stopifnot(inherits(polygon_data, 'data.frame'))
  stopifnot(inherits(covariate_data, 'data.frame'))
  stopifnot(inherits(aggregation_pixels, 'numeric'))
  stopifnot(inherits(coords, 'matrix'))
  stopifnot(inherits(startendindex, 'matrix'))
  stopifnot(inherits(mesh, 'inla.mesh'))
  
  disag_data <- list(polygon_shapefile = polygon_shapefile,
                     covariate_rasters = covariate_rasters,
                     polygon_data = polygon_data,
                     covariate_data = covariate_data,
                     aggregation_pixels = aggregation_pixels,
                     coords = coords,
                     startendindex = startendindex,
                     mesh = mesh)
  
  class(disag_data) <- c('disag.data', 'list')
  
  return(disag_data)
}
