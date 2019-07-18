context("Plotting data")

test_that("Check plot_polygon_data function works as expected", {
  
  cds1 <- rbind(c(-120,-20), c(-100,5), c(-60, 0), c(-100,-60), c(-120,-20))
  cds2 <- rbind(c(30,0), c(50,60), c(70,0), c(70,-55), c(30,0))
  cds3 <- rbind(c(10,20), c(10,50), c(-60,40), c(-30,-30), c(30,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('1', '2', '3'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  p <- plot_polygon_data(spdf)
  expect_error(plot_polygon_data(polys))
  expect_is(p, 'ggplot')
  
})

test_that("Check plot_covariate_data function works as expected", {
  
  # Create raster stack
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  cov_stack <- raster::stack(r, r)
  
  p <- plot_covariate_data(cov_stack)
  expect_error(plot_covariate_data(r))
  expect_is(p, 'trellis')
  
})

test_that("Check plot_inla_mesh function works as expected", {
  
  cds1 <- rbind(c(-120,-20), c(-100,5), c(-60, 0), c(-100,-60), c(-120,-20))
  cds2 <- rbind(c(30,0), c(50,60), c(70,0), c(70,-55), c(30,0))
  cds3 <- rbind(c(10,20), c(10,50), c(-60,40), c(-30,-30), c(30,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('1', '2', '3'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  my_mesh <- build_mesh(spdf)
  
  p <- plot_inla_mesh(my_mesh)
  expect_error(plot_inla_mesh(spdf))
  expect_is(p, 'NULL')
  
})

test_that("Check plot.disag.data function works as expected", {
  
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
  
  result <- prepare_data(polygon_shapefile = spdf, 
                         covariate_rasters = cov_rasters)
  
  p <- plot(result)
  
  expect_error(plot_polygon_data(polys))
  expect_is(p, 'list')
  expect_equal(length(p), 2)
  expect_equal(names(p), c('polygon', 'covariates'))
  
})

