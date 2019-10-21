
context("Build mesh")

test_that("build_mesh behaves as expected", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  polygons <- list()
  n_polygon_per_side <- 10
  n_polygons <- n_polygon_per_side * n_polygon_per_side
  n_pixels_per_side <- n_polygon_per_side * 2
  
  for(i in 1:n_polygons) {
    row <- ceiling(i/n_polygon_per_side)
    col <- ifelse(i %% n_polygon_per_side != 0, i %% n_polygon_per_side, n_polygon_per_side)
    xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
    polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
  }
  
  polys <- do.call(raster::spPolygons, polygons)
  
  response_df <- data.frame(area_id = 1:n_polygons, response = runif(n_polygons, min = 0, max = 1000))
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)

  my_mesh <- build_mesh(spdf)

  expect_error(build_mesh(response_df))
  expect_error(build_mesh(spdf, mesh.args = c(4, 8)))
  expect_is(my_mesh, 'inla.mesh')
  expect_is(build_mesh(spdf, mesh.args = list(max.edge = c(50, 100))), 'inla.mesh')
  
})