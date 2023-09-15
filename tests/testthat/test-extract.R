
context("Extract covariates and polygon data")

polygons <- list()
n_polygon_per_side <- 10
n_polygons <- n_polygon_per_side * n_polygon_per_side
n_pixels_per_side <- n_polygon_per_side * 2

for(i in 1:n_polygons) {
  row <- ceiling(i/n_polygon_per_side)
  col <- ifelse(i %% n_polygon_per_side != 0, i %% n_polygon_per_side, n_polygon_per_side)
  xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
  polygons[[i]] <- list(cbind(c(xmin, xmax, xmax, xmin, xmin),
                              c(ymax, ymax, ymin, ymin, ymax)))
}

polys <- lapply(polygons,sf::st_polygon)
N <- floor(runif(n_polygons, min = 1, max = 100))
response_df <- data.frame(area_id = 1:n_polygons, response = runif(n_polygons, min = 0, max = 1000))
response_binom_df <- data.frame(area_id = 1:n_polygons, response = N*runif(n_polygons, min = 0, max = 1), sample_size = N)

spdf <- sf::st_sf(response_df, geometry = polys)
spdf_binom <- sf::st_sf(response_binom_df, geometry = polys)

# Create raster stack
r <- terra::rast(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
terra::ext(r) <- terra::ext(spdf)
r[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side), 3))
r2 <- terra::rast(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
terra::ext(r2) <- terra::ext(spdf)
r2[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ceiling(x/n_pixels_per_side), 3))
cov_stack <- c(r, r2)

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
  r <- terra::rast(ncol=20, nrow=20)
  r[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
  r2 <- terra::rast(ncol=20, nrow=20)
  r2[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
  cov_stack <- c(r, r2)
  terra::writeRaster(r, paste0(tempdir(), '/cov1.tif'), overwrite = TRUE)
  terra::writeRaster(r2, paste0(tempdir(), '/cov2.tif'), overwrite = TRUE)

  expect_is(getCovariateRasters(tempdir(), '.tif$', spdf), 'SpatRaster')

})

test_that("extractCoordsForMesh function behaves as it should", {

  skip_on_cran()

  # cl <- parallel::makeCluster(2)
  # doParallel::registerDoParallel(cl)
  # cov_data <- parallelExtract(cov_stack, spdf, fun = NULL, id = )
  # parallel::stopCluster(cl)
  # foreach::registerDoSEQ()

  cov_data <- terra::extract(cov_stack, spdf, cells=TRUE, na.rm=TRUE, ID=TRUE)
  names(cov_data)[1] <- 'area_id'

  result <- extractCoordsForMesh(cov_stack, cov_data$cellid)

  result2 <- extractCoordsForMesh(cov_stack)

  expect_error(extractCoordsForMesh(cov_data$cellid, cov_stack))
  expect_is(result, 'matrix')
  expect_is(result2, 'matrix')

})
