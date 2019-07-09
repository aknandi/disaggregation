
context("Extract covariates and polygon data")

test_that("parallelExtract gives errors when it should", {
  
  cds1 <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
  cds2 <- rbind(c(80,0), c(100,60), c(120,0), c(120,-55), c(80,0))
  cds3 <- rbind(c(30,20), c(30,50), c(-40,40), c(-10,-30), c(50,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('A', 'B', 'C'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  # Create raster stack
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  r <- raster::stack(r, r)
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  
  expect_error(parallelExtract(spdf, r, fun = NULL, id = 'area_id'))
  expect_error(parallelExtract(r, spdf, fun = NULL, id = 'id'))
  
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
})

test_that("parallelExtract give the right for of output", {
  
  cds1 <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
  cds2 <- rbind(c(80,0), c(100,60), c(120,0), c(120,-55), c(80,0))
  cds3 <- rbind(c(30,20), c(30,50), c(-40,40), c(-10,-30), c(50,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)

  response_df <- data.frame(area_id = c('A', 'B', 'C'), response = c(4, 7, 2))

  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  # Create raster stack
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  r <- raster::stack(r, r)

  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  result1 <- parallelExtract(r, spdf, fun = NULL, id = 'area_id')
  result2 <- parallelExtract(r, spdf, fun = sum, id = 'area_id')
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  expect_is(result1, 'data.frame')
  expect_equal(unique(result1$area_id), spdf$area_id)
  expect_equal(ncol(result1), raster::nlayers(r) + 2)
  expect_equal(names(r), names(result1)[-c(1,2)])
  expect_equal(length(unique(result1$area_id)), length(spdf))
  
  expect_is(result2, 'data.frame')
  expect_equal(result2$area_id, spdf$area_id)
  expect_equal(ncol(result2), raster::nlayers(r) + 1)
  expect_equal(names(r), names(result2)[-c(1)])
  expect_equal(length(unique(result2$area_id)), length(spdf))
  
  # Save in tempdir() to be use in later test
  write.csv(result1, file = paste0(tempdir(), '/cov_data.csv'))
  
})

test_that("getPolygonData function", {
  
  # Create SpatialPolygonDataFrame
  cds1 <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
  cds2 <- rbind(c(80,0), c(100,60), c(120,0), c(120,-55), c(80,0))
  cds3 <- rbind(c(30,20), c(30,50), c(-40,40), c(-10,-30), c(50,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)

  response_df <- data.frame(area_id = c('A', 'B', 'C'), response = c(4, 7, 2))

  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)

  expect_error(getPolygonData(spdf, id_var = 'id', response_var = 'response'))
  expect_error(getPolygonData(spdf, id_var = 'area_id', response_var = 'data'))
  
  result <- getPolygonData(spdf, id_var = 'area_id', response_var = 'response')
  
  expect_is(result, 'data.frame')
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), nrow(spdf))
  expect_equal(result$area_id, spdf$area_id)
  expect_equal(result$response, spdf$response)
  
})

test_that("getCovariateData function gives errors when it should", {
  
  # Create SpatialPolygonDataFrame
  cds1 <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
  cds2 <- rbind(c(80,0), c(100,60), c(120,0), c(120,-55), c(80,0))
  cds3 <- rbind(c(30,20), c(30,50), c(-40,40), c(-10,-30), c(50,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('A', 'B', 'C'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  expect_error(getCovariateRasters('/home/rasters', '.tif$', spdf))
  expect_error(getCovariateRasters(tempdir(), '.tif$', spdf))
  
  # Save .tif files in tempdir()
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  raster::writeRaster(r, paste0(tempdir(), '/cov1.tif'))
  raster::writeRaster(r, paste0(tempdir(), '/cov2.tif'))
  
  expect_is(getCovariateRasters(tempdir(), '.tif$', spdf), 'RasterBrick')
  
})

test_that("extractCoordsForMesh function behaves as it should", {

  # Create raster stack
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  r <- raster::stack(r, r)
  
  cov_data <- read.csv(paste0(tempdir(), '/cov_data.csv'))

  expect_error(extractCoordsForMesh(cov_data, r))
  expect_is(extractCoordsForMesh(r, cov_data), 'matrix')

})