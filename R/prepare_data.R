#' Function to collect all data together for disaggregation model
#'
#' @param polygon_shapefile SpatialPolygonDataFrame containing the response data 
#' @param covariate_rasters RasterStack of covariates
#' @param id_var Name of column in SpatialPolygonDataFrame object with the polygon id
#' @param response_var Name of column in SpatialPolygonDataFrame object with the response data
#' @param mesh.args list of parameters that control the mesh structure with the same names as used by INLA
#'
#' @name prepare_data
#'
#' @examples 
#' \dontrun {
#'  cds1 <- rbind(c(-120,-20), c(-100,5), c(-60, 0), c(-100,-60), c(-120,-20))
#'  cds2 <- rbind(c(30,0), c(50,60), c(70,0), c(70,-55), c(30,0))
#'  cds3 <- rbind(c(10,20), c(10,50), c(-60,40), c(-30,-30), c(30,10))
#'  polys <- raster::spPolygons(cds1, cds2, cds3)
#' 
#'  response_df <- data.frame(area_id = c('1', '2', '3'), response = c(4, 7, 2))
#' 
#'  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
#' 
#'  # Create raster stack
#'  r <- raster::raster(ncol=36, nrow=18)
#'  r[] <- 1:raster::ncell(r)
#'  cov_rasters <- raster::stack(r, r)
#' 
#'  test_data <- prepare_data(polygon_shapefile = spdf, 
#'                            covariate_rasters = cov_rasters)
#' } 
#'                    
#' @export

prepare_data <- function(polygon_shapefile, 
                         covariate_rasters, 
                         id_var = 'area_id', 
                         response_var = 'response', 
                         mesh.args = NULL) {
  
  stopifnot(inherits(polygon_shapefile, 'SpatialPolygonsDataFrame'))
  stopifnot(inherits(covariate_rasters, 'Raster'))
  stopifnot(inherits(id_var, 'character'))
  stopifnot(inherits(response_var, 'character'))
  if(!is.null(mesh.args)) stopifnot(inherits(mesh.args, 'list'))
  
  polygon_data <- getPolygonData(polygon_shapefile, id_var = id_var, response_var = response_var)
  
  cl <- parallel::makeCluster(raster::nlayers(covariate_rasters))
  doParallel::registerDoParallel(cl)
  covariate_data <- parallelExtract(covariate_rasters, polygon_shapefile, fun = NULL, id = id_var)
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  coords <- extractCoordsForMesh(covariate_rasters, covariate_data)
  
  startendindex <- getStartendindex(covariate_data, polygon_data, id_var = id_var)
  
  mesh <- build_mesh(polygon_shapefile, mesh.args)
  
  disag_data <- list(polygon_data = polygon_data,
                     covariate_data = covariate_data,
                     coords = coords,
                     startendindex = startendindex,
                     mesh = mesh)
  
  class(disag_data) <- c('disag.data', 'list')
  
  return(disag_data)
  
}

#' Function to fit the disaggregation model
#'
#' @param polygon_data data.frame with two columns: polygon id and response
#' @param covariate_data data.frame with cell id, polygon id and covariate columns
#' @param coords coordinates of the covariate data points
#' @param startendindex matrix containing the start and end index for each polygon
#' @param mesh inla.mesh object to use in the fit
#' 
#' @name fit_model
#'
#' @examples 
#' \dontrun{
#'   as.disag.data(polygon_data, covariate_data, coords, startendindex, mesh)
#'  }
#' 
#' @export


as.disag.data <- function(polygon_data, covariate_data, coords, startendindex, mesh) {
  
  stopifnot(inherits(polygon_data, 'data.frame'))
  stopifnot(inherits(covariate_data, 'data.frame'))
  stopifnot(inherits(coords, 'matrix'))
  stopifnot(inherits(startendindex, 'matrix'))
  stopifnot(inherits(mesh, 'inla.mesh'))
  
  disag_data <- list(polygon_data = polygon_data,
                     covariate_data = covariate_data,
                     coords = coords,
                     startendindex = startendindex,
                     mesh = mesh)
  
  class(disag_data) <- c('disag.data', 'list')
  
  return(disag_data)
}
