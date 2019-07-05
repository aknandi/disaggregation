#' Build mesh for disaggregaton model
#'
#' @param shapes shapefile covering the region under investigation
#'
#' @importFrom maptools unionSpatialPolygons
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

build_mesh <- function(shapes) {

  outline <- unionSpatialPolygons(shapes, IDs = rep(1, length(shapes)))

  coords <- list()
  for(i in seq_len(length(outline@polygons[[1]]@Polygons))){
    coords[[i]] <- outline@polygons[[1]]@Polygons[[i]]@coords
  }
  coords <- do.call(rbind, coords)

  outline.hull <- INLA::inla.nonconvex.hull(coords, convex = -0.01, concave = -0.5, resolution = 300)
  mesh <- INLA::inla.mesh.2d(boundary = outline.hull, max.edge = c(3, 8), cutoff = 0.4, offset = c(1, 15))

  return(mesh)
}
