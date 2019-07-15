
context("Extract covariates and polygon data")

test_that("parallelExtract gives errors when it should", {
  
  spdf <- rgdal::readOGR(paste0(tempdir(), '/test_spdf.shp'))
  
  # Create raster stack
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  r <- raster::stack(r, r)
  raster::writeRaster(r, paste0(tempdir(), '/test_cov_stack.tif'))
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  
  expect_error(parallelExtract(spdf, r, fun = NULL, id = 'area_id'))
  expect_error(parallelExtract(r, spdf, fun = NULL, id = 'id'))
  
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
})

test_that("parallelExtract give the right for of output", {
  
  spdf <- rgdal::readOGR(paste0(tempdir(), '/test_spdf.shp'))
  
  r <- raster::stack(paste0(tempdir(), '/test_cov_stack.tif'))

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
  write.csv(result1, file = paste0(tempdir(), '/test_cov_data.csv'))
  
})

test_that("getPolygonData function", {
  
  spdf <- rgdal::readOGR(paste0(tempdir(), '/test_spdf.shp'))

  expect_error(getPolygonData(spdf, id_var = 'id', response_var = 'response'))
  expect_error(getPolygonData(spdf, id_var = 'area_id', response_var = 'data'))
  
  result <- getPolygonData(spdf, id_var = 'area_id', response_var = 'response')
  
  write.csv(result, file = paste0(tempdir(), '/test_polygon_data.csv'))
  
  expect_is(result, 'data.frame')
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), nrow(spdf))
  expect_equal(result$area_id, spdf$area_id)
  expect_equal(result$response, spdf$response)
  
})

test_that("getCovariateData function gives errors when it should", {
  
  spdf <- rgdal::readOGR(paste0(tempdir(), '/test_spdf.shp'))
  
  expect_error(getCovariateRasters('/home/rasters', '.tif$', spdf))
  
  # Save .tif files in tempdir()
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  raster::writeRaster(r, paste0(tempdir(), '/cov1.tif'))
  raster::writeRaster(r, paste0(tempdir(), '/cov2.tif'))
  
  expect_is(getCovariateRasters(tempdir(), '.tif$', spdf), 'RasterBrick')
  
})

test_that("extractCoordsForMesh function behaves as it should", {

  r <- raster::stack(paste0(tempdir(), '/test_cov_stack.tif'))
  cov_data <- read.csv(paste0(tempdir(), '/test_cov_data.csv'))

  result <- extractCoordsForMesh(r, cov_data)
  save(result, file = paste0(tempdir(), '/test_coords.RData'))
  
  expect_error(extractCoordsForMesh(cov_data, r))
  expect_is(result, 'matrix')

})