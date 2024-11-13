context("Preparing data")

test_that("Check prepare_data function works as expected", {

  result <- prepare_data(polygon_shapefile = spdf,
                         covariate_rasters = cov_stack)

  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data',
                                'aggregation_pixels', 'coords_for_fit', 'coords_for_prediction', 'start_end_index', 'mesh'))
  expect_is(result$polygon_shapefile, 'sf')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, 'SpatRaster')
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coords_for_fit, 'matrix')
  expect_is(result$coords_for_prediction, 'matrix')
  expect_is(result$start_end_index, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(sum(is.na(result$polygon_data$N)), length(result$polygon_data$N))
  expect_equal(nrow(result$polygon_data), nrow(result$start_end_index))
  expect_equal(nrow(result$covariate_data), nrow(result$coords_for_fit))

})

test_that("Check prepare_data function with sample size works as expected", {

  result <- prepare_data(polygon_shapefile = spdf_binom,
                         covariate_rasters = cov_stack,
                         sample_size_var = 'sample_size',
                         make_mesh = FALSE)

  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data',
                                'aggregation_pixels', 'coords_for_fit', 'coords_for_prediction', 'start_end_index', 'mesh'))
  expect_is(result$polygon_shapefile, 'sf')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, 'SpatRaster')
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coords_for_fit, 'matrix')
  expect_is(result$coords_for_prediction, 'matrix')
  expect_is(result$start_end_index, 'matrix')
  expect_true(is.null(result$mesh))
  expect_equal(sum(is.na(result$polygon_data$N)), 0)
  expect_equal(nrow(result$polygon_data), nrow(result$start_end_index))
  expect_equal(nrow(result$covariate_data), nrow(result$coords_for_fit))

})

test_that("Check prepare_data function deals with NAs as expected", {

  cov_stack_na <- cov_stack
  cov_stack_na[[1]][c(1:10)] <- NA

  aggregation_raster_na <- r
  aggregation_raster_na[c(1:10)] <- NA

  aggregation_raster_zero <- r
  aggregation_raster_zero[c(1:2, 21:22)] <- 0

  expect_error(prepare_data(polygon_shapefile = spdf_na, covariate_rasters = cov_stack, make_mesh = FALSE))
  expect_error(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_stack_na, make_mesh = FALSE))
  expect_error(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_stack, aggregation_raster = aggregation_raster_na, make_mesh = FALSE))
  expect_error(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_stack, aggregation_raster = aggregation_raster_zero, make_mesh = FALSE))

  result <- prepare_data(polygon_shapefile = spdf_na,
                         covariate_rasters = cov_stack_na,
                         aggregation_raster = aggregation_raster_na,
                         na_action = TRUE,
                         make_mesh = FALSE)

  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data',
                                'aggregation_pixels', 'coords_for_fit', 'coords_for_prediction', 'start_end_index', 'mesh'))
  expect_is(result$polygon_shapefile, 'sf')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, 'SpatRaster')
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coords_for_fit, 'matrix')
  expect_is(result$start_end_index, 'matrix')
  expect_true(is.null(result$mesh))
  expect_equal(nrow(result$polygon_data), nrow(result$start_end_index))
  expect_equal(nrow(result$covariate_data), nrow(result$coords_for_fit))
  expect_equal(sum(is.na(result$polygon_data$response)), 0)
  expect_equal(sum(is.na(result$covariate_data)), 0)
  expect_equal(sum(is.na(result$aggregation_pixels)), 0)
  expect_equal(nrow(result$polygon_shapefile), nrow(spdf_na) - 1)

  result <- prepare_data(polygon_shapefile = spdf_na,
                         covariate_rasters = cov_stack_na,
                         aggregation_raster = aggregation_raster_zero,
                         na_action = TRUE,
                         make_mesh = FALSE)

  expect_equal(nrow(result$polygon_shapefile), nrow(spdf_na) - 2)
})


test_that("Check as.disag_data function works as expected", {

  polygon_data <- getPolygonData(spdf, id_var = 'area_id', response_var = 'response')

  cov_data <- terra::extract(cov_stack, spdf, cells=TRUE, na.rm=TRUE, ID=TRUE)
  names(cov_data)[1] <- 'area_id'

  aggregation_data <- rep(1, nrow(cov_data))

  coords_for_fit <- extractCoordsForMesh(cov_stack, cov_data$cellid)

  coords_for_prediction <- extractCoordsForMesh(cov_stack)

  start_end_index <- getStartendindex(cov_data, polygon_data, 'area_id')

  result <- as.disag_data(spdf,
                          list('area_id', 'response'),
                          cov_stack,
                          polygon_data,
                          cov_data,
                          aggregation_data,
                          coords_for_fit,
                          coords_for_prediction,
                          start_end_index,
                          mesh = NULL)

  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data',
                                'aggregation_pixels', 'coords_for_fit', 'coords_for_prediction', 'start_end_index', 'mesh'))
  expect_is(result$polygon_shapefile, 'sf')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, 'SpatRaster')
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coords_for_fit, 'matrix')
  expect_is(result$coords_for_prediction, 'matrix')
  expect_is(result$start_end_index, 'matrix')
  expect_true(is.null(result$mesh))
  expect_equal(nrow(result$polygon_data), nrow(result$start_end_index))
  expect_equal(nrow(result$covariate_data), nrow(result$coords_for_fit))

})

test_that("Check prepare_data warns about non-numeric covariates", {
  r3 <- terra::rast(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
  terra::ext(r3) <- terra::ext(spdf)
  r3[] <- sapply(1:terra::ncell(r), function(x) as.integer(x))

  cov_stack <- c(r, r2, r3)
  names(cov_stack) <- c('layer1', 'layer2', 'layer3')

  expect_warning(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_stack),
                 "The values of layer3 are not numeric")

})
