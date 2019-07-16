
context("Build mesh")

test_that("build_mesh behaves as expected", {
  
  cds1 <- rbind(c(-120,-20), c(-100,5), c(-60, 0), c(-100,-60), c(-120,-20))
  cds2 <- rbind(c(30,0), c(50,60), c(70,0), c(70,-55), c(30,0))
  cds3 <- rbind(c(10,20), c(10,50), c(-60,40), c(-30,-30), c(30,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('1', '2', '3'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)

  my_mesh <- build_mesh(spdf)

  expect_error(build_mesh(response_df))
  expect_error(build_mesh(spdf, mesh.args = c(4, 8)))
  expect_is(my_mesh, 'inla.mesh')
  expect_is(build_mesh(spdf, mesh.args = list(max.edge = c(50, 100))), 'inla.mesh')
  
})