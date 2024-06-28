
context("Extract covariates and polygon data")

test_that("getPolygonData function", {

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

  expect_error(getCovariateRasters('/home/rasters', '.tif$', spdf))

  # Save .tif files in tempdir()
  terra::writeRaster(r, paste0(tempdir(), '/cov1.tif'), overwrite = TRUE)
  terra::writeRaster(r2, paste0(tempdir(), '/cov2.tif'), overwrite = TRUE)

  expect_is(getCovariateRasters(tempdir(), '.tif$', spdf), 'SpatRaster')

})

test_that("extractCoordsForMesh function behaves as it should", {

  cov_data <- terra::extract(cov_stack, spdf, cells = TRUE, na.rm = TRUE, ID = TRUE)
  names(cov_data)[1] <- 'area_id'

  result <- extractCoordsForMesh(cov_stack, cov_data$cell)

  result2 <- extractCoordsForMesh(cov_stack)

  expect_error(extractCoordsForMesh(cov_data$cellid, cov_stack))
  expect_is(result, 'matrix')
  expect_is(result2, 'matrix')

})
