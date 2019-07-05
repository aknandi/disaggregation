
context("Build mesh")

test_that("build_mesh behaves as expected", {
  
  cds1 <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
  cds2 <- rbind(c(80,0), c(100,60), c(120,0), c(120,-55), c(80,0))
  cds3 <- rbind(c(30,20), c(30,50), c(-40,40), c(-10,-30), c(50,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('A', 'B', 'C'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  expect_error(build_mesh(response_df))
  expect_is(build_mesh(spdf), 'inla.mesh')
  
})