#' Build mesh for disaggregaton model
#'
#' @param shapes shapefile covering the region under investigation
#' @param mesh.args list of parameters that control the mesh structure with the same names as used by INLA
#'
#' @name build_mesh
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
#'  my_mesh <- build_mesh(spdf)
#' }
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
