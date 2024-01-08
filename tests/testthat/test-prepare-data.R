context("Preparing data")

test_that("Check prepare_data function works as expected", {

  skip_on_cran()

  result <- prepare_data(polygon_shapefile = spdf,
                         covariate_rasters = cov_stack)

  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data',
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'sf')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, 'SpatRaster')
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$coordsForPrediction, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(sum(is.na(result$polygon_data$N)), length(result$polygon_data$N))
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))

})

test_that("Check prepare_data function with sample size works as expected", {

  skip_on_cran()

  result <- prepare_data(polygon_shapefile = spdf_binom,
                         covariate_rasters = cov_stack,
                         sample_size_var = 'sample_size',
                         makeMesh = FALSE)

  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data',
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'sf')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, 'SpatRaster')
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$coordsForPrediction, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_true(is.null(result$mesh))
  expect_equal(sum(is.na(result$polygon_data$N)), 0)
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))

})

test_that("Check prepare_data function deals with NAs as expected", {

  skip_on_cran()

  cov_stack_na <- cov_stack
  cov_stack_na[[1]][c(1:10)] <- NA

  aggregation_raster_na <- r
  aggregation_raster_na[c(1:10)] <- NA

  expect_error(prepare_data(polygon_shapefile = spdf_na, covariate_rasters = cov_stack, makeMesh = FALSE))
  expect_error(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_stack_na, makeMesh = FALSE))
  expect_error(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_stack, aggregation_raster = aggregation_raster_na, makeMesh = FALSE))

  result <- prepare_data(polygon_shapefile = spdf_na,
                         covariate_rasters = cov_stack_na,
                         aggregation_raster = aggregation_raster_na,
                         na.action = TRUE,
                         makeMesh = FALSE)

  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data',
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'sf')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, 'SpatRaster')
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_true(is.null(result$mesh))
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))
  expect_equal(sum(is.na(result$polygon_data$response)), 0)
  expect_equal(sum(is.na(result$covariate_data)), 0)
  expect_equal(sum(is.na(result$aggregation_pixels)), 0)
  expect_equal(nrow(result$polygon_shapefile), nrow(spdf_na) - 1)
})


test_that("Check as.disag_data function works as expected", {

  skip_on_cran()

  polygon_data <- getPolygonData(spdf, id_var = 'area_id', response_var = 'response')

  cov_data <- terra::extract(cov_stack, spdf, cells=TRUE, na.rm=TRUE, ID=TRUE)
  names(cov_data)[1] <- 'area_id'

  aggregation_data <- rep(1, nrow(cov_data))

  coordsForFit <- extractCoordsForMesh(cov_stack, cov_data$cellid)

  coordsForPrediction <- extractCoordsForMesh(cov_stack)

  startendindex <- getStartendindex(cov_data, polygon_data, 'area_id')

  result <- as.disag_data(spdf,
                          list('area_id', 'response'),
                          cov_stack,
                          polygon_data,
                          cov_data,
                          aggregation_data,
                          coordsForFit,
                          coordsForPrediction,
                          startendindex,
                          mesh = NULL)

  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data',
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'sf')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, 'SpatRaster')
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$coordsForPrediction, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_true(is.null(result$mesh))
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))

})

