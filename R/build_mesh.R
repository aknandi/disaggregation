#' Build mesh for disaggregaton model
#' 
#' \emph{build_mesh} function takes a SpatialPolygons object and mesh arguments to build an appropriate mesh for the spatial field.
#'
#' The mesh is created by finding a tight boundary around the polygon data, and creating a fine mesh within the boundary 
#' and a coarser mesh outside. This speeds up computation time by only having a very fine mesh within the area of interest 
#' and having a small region outside with a coarser mesh to avoid edge effects.
#' 
#' Six mesh parameters can be specified as arguments: \emph{convex}, \emph{concave} and \emph{resolution},
#' to control the boundary of the inner mesh, and \emph{max.edge}, \emph{cut} and \emph{offset}, to control the  mesh itself,
#' with the names meaning the same as used by INLA functions \emph{inla.convex.hull} and \emph{inla.mesh.2d}.
#' 
#' Defaults are:
#' pars <- list(convex = -0.01, concave = -0.5, resolution = 300, max.edge = c(3.0, 8),  cut = 0.4, offset = c(1, 15)).
#' 
#' @param shapes shapefile covering the region under investigation.
#' @param mesh.args list of parameters that control the mesh structure. \emph{convex}, \emph{concave} and \emph{resolution},
#' to control the boundary of the inner mesh, and \emph{max.edge}, \emph{cut} and \emph{offset}, to control the  mesh itself,
#' with the parameters having the same meaning as in the INLA functions \emph{inla.convex.hull} and \emph{inla.mesh.2d}.
#'
#' @return An inla.mesh object
#' 
#' @name build_mesh
#'
#' @examples 
#' \dontrun{
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
#'  my_mesh <- build_mesh(spdf)
#' }
#'
#'
#'
#' @export

build_mesh <- function(shapes, mesh.args = NULL) {

  stopifnot(inherits(shapes, 'SpatialPolygons'))
  if(!is.null(mesh.args)) stopifnot(inherits(mesh.args, 'list'))
  
  limits <- sp::bbox(shapes)
  hypotenuse <- sqrt((limits[1,2] - limits[1,1])^2 + (limits[2,2] - limits[2,1])^2)
  maxedge <- hypotenuse/10
  
  
  pars <- list(convex = -0.01,
               concave = -0.5,
               resolution = 300,
               max.edge = c(maxedge, maxedge * 2), 
               cut = 0.1, 
               offset = c(hypotenuse / 10, hypotenuse / 10))

  
  pars[names(mesh.args)] <- mesh.args

  outline <- sf::st_union(sf::st_as_sf(shapes))
  coords <- sf::st_coordinates(outline)
  
  

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
