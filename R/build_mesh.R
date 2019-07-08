#' Build mesh for disaggregaton model
#'
#' @param shapes shapefile covering the region under investigation
#' @param mesh.args list of parameters that control the mesh structure with the same names as used by INLA
#'
#' @name build_mesh
#'
#' @examples
#'  \dontrun{
#'   build_mesh(my_shapefile)
#'  }
#'
#'
#' @export

build_mesh <- function(shapes, mesh.args = NULL) {

  stopifnot(inherits(shapes, 'SpatialPolygons'))
  if(!is.null(mesh.args)) stopifnot(inherits(mesh.args, 'list'))
  
  pars <- list(convex = -0.01,
               concave = -0.5,
               resolution = 300,
               max.edge = c(3.0, 8), 
               cut = 0.4, 
               offset = c(1, 15))
  
  pars[names(mesh.args)] <- mesh.args

  outline <- maptools::unionSpatialPolygons(shapes, IDs = rep(1, length(shapes)))

  coords <- list()
  for(i in seq_len(length(outline@polygons[[1]]@Polygons))){
    coords[[i]] <- outline@polygons[[1]]@Polygons[[i]]@coords
  }
  coords <- do.call(rbind, coords)

  outline.hull <- INLA::inla.nonconvex.hull(coords, 
                                            convex = pars$convex, 
                                            concave = pars$concave,
                                            resolution = pars$resolution)
  
  mesh <- INLA::inla.mesh.2d( 
    boundary = outline.hull,
    max.edge = pars$max.edge, 
    cut = pars$cut, 
    offset = pars$offset)
  
  
  return(mesh)
}
