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
#'  \dontrun{
#'   prepare_data(my_shapefile, cov_path)
#'  }
#'
#' @export

prepare_data <- function(polygon_shapefile, 
                         covariate_rasters, 
                         id_var = 'area_id', 
                         response_var = 'response', 
                         mesh.args = NULL) {
  
  polygon_data <- getPolygonData(polygon_shapefile, id_var = id_var, response_var = response_var)
  
  covariate_data <- parallelExtract(covariate_rasters, polygon_shapefile, fun = NULL, id = id_var)
  
  coords <- extractCoordsForMesh(covariate_rasters, covariate_data)
  
  startendindex <- getStartendindex(covariate_data, polygon_data, id_var = id_var)
  
  mesh <- build_mesh(polygon_shapefile, mesh.args)
  
  return(list(polygon_data = polygon_data,
              covariate_data = covariate_data,
              coords = coords,
              startendindex = startendindex,
              mesh = mesh))
  
}
