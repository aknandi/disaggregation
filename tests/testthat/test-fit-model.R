
context("Fitting model")

test_that("fit_model produces errors whe expected", {
  
  cds1 <- rbind(c(-120,-20), c(-100,5), c(-60, 0), c(-100,-60), c(-120,-20))
  cds2 <- rbind(c(30,0), c(50,60), c(70,0), c(70,-55), c(30,0))
  cds3 <- rbind(c(10,20), c(10,50), c(-60,40), c(-30,-30), c(30,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('1', '2', '3'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  # Create raster stack
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  cov_rasters <- raster::stack(r, r)
  
  test_data <- prepare_data(polygon_shapefile = spdf, 
                            covariate_rasters = cov_rasters,
                            mesh.args = list(max.edge = c(50, 100)))
  
  result <- fit_model(test_data, its = 2)
  
  expect_error(fit_model(list()))
  expect_error(fit_model(test_data, its = 'its'))
  expect_is(result, 'fit.result')
  expect_equal(length(result), 3)
  
})