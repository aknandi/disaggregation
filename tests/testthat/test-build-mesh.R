
context("Build mesh")

test_that("build_mesh behaves as expected", {
  
  polygons <- list()
  for(i in 1:100) {
    row <- ceiling(i/10)
    col <- ifelse(i %% 10 != 0, i %% 10, 10)
    xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
    polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
  }
  
  polys <- do.call(raster::spPolygons, polygons)
  response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 1000))
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)

  my_mesh <- build_mesh(spdf)

  expect_error(build_mesh(response_df))
  expect_error(build_mesh(spdf, mesh.args = c(4, 8)))
  expect_is(my_mesh, 'inla.mesh')
  expect_is(build_mesh(spdf, mesh.args = list(max.edge = c(50, 100))), 'inla.mesh')
  
})