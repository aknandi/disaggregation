context("Predict model")

test_that("Check predict.disag_model function works as expected", {

  pred2 <- predict(result)

  expect_is(pred2, 'disag_prediction')
  expect_equal(length(pred2), 2)
  expect_equal(names(pred2), c('mean_prediction', 'uncertainty_prediction'))

  expect_is(pred2$mean_prediction, 'list')
  expect_equal(length(pred2$mean_prediction), 4)
  expect_is(pred2$mean_prediction$prediction, 'SpatRaster')
  expect_is(pred2$mean_prediction$field, 'SpatRaster')
  expect_true(is.null(pred2$mean_prediction$iid))
  expect_is(pred2$mean_prediction$covariates, 'SpatRaster')

  expect_is(pred2$uncertainty_prediction, 'list')
  expect_equal(length(pred2$uncertainty_prediction), 2)
  expect_equal(names(pred2$uncertainty_prediction), c('realisations', 'predictions_ci'))
  expect_is(pred2$uncertainty_prediction$realisations, 'SpatRaster')
  expect_is(pred2$uncertainty_prediction$predictions_ci, 'SpatRaster')
  expect_equal(terra::nlyr(pred2$uncertainty_prediction$realisations), 100)
  expect_equal(terra::nlyr(pred2$uncertainty_prediction$predictions_ci), 2)

  pred2 <- predict(result, predict_iid = TRUE, N = 10)

  expect_is(pred2, 'disag_prediction')
  expect_equal(length(pred2), 2)
  expect_equal(names(pred2), c('mean_prediction', 'uncertainty_prediction'))

  expect_is(pred2$mean_prediction, 'list')
  expect_equal(length(pred2$mean_prediction), 4)
  expect_equal(names(pred2$mean_prediction), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(pred2$mean_prediction$prediction, 'SpatRaster')
  expect_is(pred2$mean_prediction$field, 'SpatRaster')
  expect_is(pred2$mean_prediction$iid, 'SpatRaster')
  expect_is(pred2$mean_prediction$covariates, 'SpatRaster')

  expect_is(pred2$uncertainty_prediction, 'list')
  expect_equal(length(pred2$uncertainty_prediction), 2)
  expect_equal(names(pred2$uncertainty_prediction), c('realisations', 'predictions_ci'))
  expect_is(pred2$uncertainty_prediction$realisations, 'SpatRaster')
  expect_is(pred2$uncertainty_prediction$predictions_ci, 'SpatRaster')
  expect_equal(terra::nlyr(pred2$uncertainty_prediction$realisations), 10)
  expect_equal(terra::nlyr(pred2$uncertainty_prediction$predictions_ci), 2)


  # For a model with no field or iid

  result2 <- disag_model(test_data, iterations = 100, field = FALSE, iid = FALSE)

  pred2 <- predict(result2)

  expect_is(pred2, 'disag_prediction')
  expect_equal(length(pred2), 2)
  expect_equal(names(pred2), c('mean_prediction', 'uncertainty_prediction'))

  expect_is(pred2$mean_prediction, 'list')
  expect_equal(length(pred2$mean_prediction), 4)
  expect_is(pred2$mean_prediction$prediction, 'SpatRaster')
  expect_true(is.null(pred2$mean_prediction$field))
  expect_true(is.null(pred2$mean_prediction$iid))
  expect_is(pred2$mean_prediction$covariates, 'SpatRaster')

  expect_is(pred2$uncertainty_prediction, 'list')
  expect_equal(length(pred2$uncertainty_prediction), 2)
  expect_equal(names(pred2$uncertainty_prediction), c('realisations', 'predictions_ci'))
  expect_is(pred2$uncertainty_prediction$realisations, 'SpatRaster')
  expect_is(pred2$uncertainty_prediction$predictions_ci, 'SpatRaster')
  expect_equal(terra::nlyr(pred2$uncertainty_prediction$realisations), 100)
  expect_equal(terra::nlyr(pred2$uncertainty_prediction$predictions_ci), 2)

})



test_that("Check predict.disag_model function works with new data", {

  new_data <- terra::crop(c(r, r2), c(0, 10, 0, 10))
  names(new_data) <- c('layer1', 'layer2')
  pred1 <- predict(result)
  pred2 <- predict(result, new_data, predict_iid = TRUE, N = 5)

  expect_is(pred2, 'disag_prediction')
  expect_equal(length(pred2), 2)
  expect_equal(names(pred2), c('mean_prediction', 'uncertainty_prediction'))

  expect_is(pred2$mean_prediction, 'list')
  expect_equal(length(pred2$mean_prediction), 4)
  expect_equal(names(pred2$mean_prediction), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(pred2$mean_prediction$prediction, 'SpatRaster')
  expect_true(!is.null(pred2$mean_prediction$field))
  expect_is(pred2$mean_prediction$iid, 'SpatRaster')
  expect_is(pred2$mean_prediction$covariates, 'SpatRaster')

  expect_is(pred2$uncertainty_prediction, 'list')
  expect_equal(length(pred2$uncertainty_prediction), 2)
  expect_equal(names(pred2$uncertainty_prediction), c('realisations', 'predictions_ci'))
  expect_is(pred2$uncertainty_prediction$realisations, 'SpatRaster')
  expect_is(pred2$uncertainty_prediction$predictions_ci, 'SpatRaster')
  expect_equal(terra::nlyr(pred2$uncertainty_prediction$realisations), 5)
  expect_equal(terra::nlyr(pred2$uncertainty_prediction$predictions_ci), 2)

  expect_false(identical(terra::ext(pred1$mean_prediction$prediction), terra::ext(pred2$mean_prediction$prediction)))
  expect_false(identical(terra::ext(pred1$uncertainty_prediction$realisations), terra::ext(pred2$uncertainty_prediction$realisations)))

})

test_that('Check that check_new_data works', {

  new_data <- terra::crop(c(r, r2), c(0, 10, 0, 10))
  names(new_data) <- c('layer1', 'layer2')

  nd1 <- check_new_data(new_data, result)
  expect_is(nd1, 'SpatRaster')

  nn <- new_data[[1]]
  names(nn) <- 'extra_unneeded'
  new_data2 <- c(new_data, nn)
  expect_error(check_new_data(new_data2, result), NA)

  new_data3 <- new_data[[1]]
  expect_error(check_new_data(new_data3, result), 'All covariates')

  new_data4 <- result$data$covariate_data
  expect_error(check_new_data(new_data4, result), 'new_data should be NULL or')


})

test_that('Check that setup_objects works', {

  objects <- setup_objects(result)

  expect_is(objects, 'list')
  expect_equal(length(objects), 3)
  expect_equal(names(objects), c('covariates', 'field_objects', 'iid_objects'))
  expect_is(objects$field_objects, 'list')
  expect_true(is.null(objects$iid_objects))

  new_data <- terra::crop(c(r, r2), c(0, 180, -90, 90))
  names(new_data) <- c('layer1', 'layer2')
  objects2 <- setup_objects(result, new_data)

  expect_is(objects2, 'list')
  expect_equal(length(objects2), 3)
  expect_equal(names(objects2), c('covariates', 'field_objects', 'iid_objects'))
  expect_is(objects2$field_objects, 'list')
  expect_true(is.null(objects$iid_objects))

  objects3 <- setup_objects(result, predict_iid = TRUE)

  expect_is(objects3, 'list')
  expect_equal(length(objects3), 3)
  expect_equal(names(objects3), c('covariates', 'field_objects', 'iid_objects'))
  expect_is(objects3$field_objects, 'list')
  expect_is(objects3$iid_objects, 'list')

})

test_that('Check that predict_single_raster works', {

  objects <- setup_objects(result)

  pars <- result$obj$env$last.par.best
  pars <- split(pars, names(pars))

  pred2 <- predict_single_raster(pars,
                                objects = objects,
                                link_function = result$model_setup$link)

  expect_is(pred2, 'list')
  expect_equal(length(pred2), 4)
  expect_equal(names(pred2), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(pred2$prediction, 'SpatRaster')
  expect_is(pred2$field, 'SpatRaster')
  expect_true(is.null(pred2$iid))
  expect_is(pred2$covariates, 'SpatRaster')

  objects2 <- setup_objects(result, predict_iid = TRUE)

  pred2 <- predict_single_raster(pars,
                                objects = objects2,
                                link_function = result$model_setup$link)

  expect_is(pred2, 'list')
  expect_equal(length(pred2), 4)
  expect_equal(names(pred2), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(pred2$prediction, 'SpatRaster')
  expect_is(pred2$field, 'SpatRaster')
  expect_is(pred2$iid, 'SpatRaster')
  expect_is(pred2$covariates, 'SpatRaster')

})


