
context("Preparing data")

test_that("Check prepare_data function works as expected", {
  
  # Create SpatialPolygonDataFrame
  cds1 <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
  cds2 <- rbind(c(80,0), c(100,60), c(120,0), c(120,-55), c(80,0))
  cds3 <- rbind(c(30,20), c(30,50), c(-40,40), c(-10,-30), c(50,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('A', 'B', 'C'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  # Create raster stack
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  cov_rasters <- raster::stack(r, r)
  
  result <- prepare_data(polygon_shapefile = spdf, 
                         covariate_rasters = cov_rasters, 
                         mesh.args = list(max.edge = c(4.0, 8.0)))
  
  expect_is(result, 'list')
  expect_equal(length(result), 5)
  expect_equal(names(result), c('polygon_data', 'covariate_data', 'coords', 'startendindex', 'mesh'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$coords, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coords))
  
})

