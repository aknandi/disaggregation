context("Summary functions")

test_that("Check summary.disag_data function works as expected", {

  data_summary <- summary(test_data)

  expect_is(data_summary, 'list')
  expect_equal(length(data_summary), 3)
  expect_equal(names(data_summary), c('number_polygons', 'number_covariates', 'covariate_summary'))
  expect_is(data_summary$number_polygons, 'integer')
  expect_is(data_summary$number_covariates, 'integer')
  expect_is(data_summary$covariate_summary, 'table')
  expect_equal(ncol(data_summary$covariate_summary), data_summary$number_covariates)

})

test_that("Check print.disag_data function works as expected", {

  print_output <- print(test_data)

  expect_is(print_output, 'disag_data')
  expect_equal(print_output, test_data)

})

test_that("Check summary.disag_model function works as expected", {

  model_summary <- summary(result)

  expect_is(model_summary, 'list')
  expect_equal(length(model_summary), 3)
  expect_equal(names(model_summary), c('model_params', 'nll', 'metrics'))
  expect_true(c("layer1" %in% rownames(model_summary$model_params)))
  expect_true(c("layer2" %in% rownames(model_summary$model_params)))
  expect_is(model_summary$model_params, 'matrix')
  expect_is(model_summary$nll, 'numeric')
  expect_is(model_summary$metrics, 'data.frame')
  expect_equal(dim(model_summary$metrics), c(1, 5))

})

test_that("Check print.disag_model function works as expected", {

  print_output <- print(result)

  expect_is(print_output, 'disag_model')
  expect_equal(print_output, result)

})

test_that("Check summary.disag_predictions function works as expected", {

  pred <- predict(result)

  model_summary <- summary(pred)

  expect_is(model_summary, 'list')
  expect_equal(length(model_summary), 3)
  expect_equal(names(model_summary), c('number_realisations', 'range_mean_values', 'range_iqr_values'))
  expect_is(model_summary$number_realisations, 'integer')
  expect_is(model_summary$range_mean_values, 'numeric')
  expect_is(model_summary$range_iqr_values, 'numeric')

})

test_that("Check print.disag_predictions function works as expected", {

  pred <- predict(result)

  print_output <- print(pred)

  expect_is(print_output, 'disag_prediction')
  expect_equal(print_output, pred)

})
