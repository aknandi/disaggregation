
context("Extract covariates and polygon data")

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
N <- floor(runif(n_polygons, min = 1, max = 100))
response_df <- data.frame(area_id = 1:n_polygons, response = runif(n_polygons, min = 0, max = 1000))
response_binom_df <- data.frame(area_id = 1:n_polygons, response = N*runif(n_polygons, min = 0, max = 1), sample_size = N)

spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
spdf_binom <- sp::SpatialPolygonsDataFrame(polys, response_binom_df)

# Create raster stack
r <- raster::raster(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
r <- raster::setExtent(r, raster::extent(spdf))
r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side), 3))
r2 <- raster::raster(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
r2 <- raster::setExtent(r2, raster::extent(spdf))
r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/n_pixels_per_side), 3))
cov_stack <- raster::stack(r, r2)

test_that("parallelExtract gives errors when it should", {
  
  skip_on_cran()
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  
  expect_error(parallelExtract(spdf, cov_stack, fun = NULL, id = 'area_id'))
  expect_error(parallelExtract(cov_stack, spdf, fun = NULL, id = 'id'))
  
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
})

test_that("parallelExtract give the right form of output", {
  
  skip_on_cran()
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  cov_data <- parallelExtract(cov_stack, spdf, fun = NULL, id = 'area_id')
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  expect_is(cov_data, 'data.frame')
  expect_equal(sort(as.numeric(unique(cov_data$area_id))), spdf$area_id)#
  expect_equal(ncol(cov_data), raster::nlayers(cov_stack) + 2)#
  expect_equal(names(cov_stack), names(cov_data)[-c(1,2)])#
  expect_equal(length(unique(cov_data$area_id)), length(spdf))
  
})

test_that("getPolygonData function", {
  
  skip_on_cran()
  
  expect_error(getPolygonData(spdf, id_var = 'id', response_var = 'response'))
  expect_error(getPolygonData(spdf, id_var = 'area_id', response_var = 'data'))
  
  result <- getPolygonData(spdf, id_var = 'area_id', response_var = 'response')
  result_binom <- getPolygonData(spdf_binom, id_var = 'area_id', response_var = 'response', sample_size_var = 'sample_size')
  
  expect_is(result, 'data.frame')
  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), nrow(spdf))
  expect_equal(result$area_id, spdf$area_id)
  expect_equal(result$response, spdf$response)
  expect_equal(result$N, rep(NA, nrow(result)))
  
  expect_is(result_binom, 'data.frame')
  expect_equal(ncol(result_binom), 3)
  expect_equal(nrow(result_binom), nrow(spdf_binom))
  expect_equal(result_binom$area_id, spdf_binom$area_id)
  expect_equal(result_binom$response, spdf_binom$response)
  expect_equal(result_binom$N, spdf_binom$sample_size)
  
})

test_that("getCovariateData function gives errors when it should", {
  
  skip_on_cran()
  
  expect_error(getCovariateRasters('/home/rasters', '.tif$', spdf))
  
  # Save .tif files in tempdir()
  r <- raster::raster(ncol=20, nrow=20)
  r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
  r2 <- raster::raster(ncol=20, nrow=20)
  r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
  cov_stack <- raster::stack(r, r2)
  raster::writeRaster(r, paste0(tempdir(), '/cov1.tif'), overwrite = TRUE)
  raster::writeRaster(r2, paste0(tempdir(), '/cov2.tif'), overwrite = TRUE)
  
  expect_is(getCovariateRasters(tempdir(), '.tif$', spdf), 'RasterBrick')
  
})

test_that("extractCoordsForMesh function behaves as it should", {

  skip_on_cran()
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  cov_data <- parallelExtract(cov_stack, spdf, fun = NULL, id = 'area_id')
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  result <- extractCoordsForMesh(cov_stack, cov_data$cellid)
  
  result2 <- extractCoordsForMesh(cov_stack)

  expect_error(extractCoordsForMesh(cov_data$cellid, cov_stack))
  expect_is(result, 'matrix')
  expect_is(result2, 'matrix')

})
