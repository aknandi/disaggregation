context("Summary functions")

test_that("Check summary.disag_data function works as expected", {

  skip_on_cran()

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

  skip_on_cran()

  print_output <- print(test_data)

  expect_is(print_output, 'disag_data')
  expect_equal(print_output, test_data)

})

test_that("Check summary.disag_model function works as expected", {

  skip_if_not_installed('INLA')
  skip_on_cran()

  result <- disag_model(test_data, field = FALSE, iterations = 2)

  model_summary <- summary(result)

  expect_is(model_summary, 'list')
  expect_equal(length(model_summary), 3)
  expect_equal(names(model_summary), c('model_params', 'nll', 'metrics'))
  expect_is(model_summary$model_params, 'matrix')
  expect_is(model_summary$nll, 'numeric')
  expect_is(model_summary$metrics, 'data.frame')
  expect_equal(dim(model_summary$metrics), c(1, 5))

})

test_that("Check print.disag_model function works as expected", {

  skip_if_not_installed('INLA')
  skip_on_cran()

  result <- disag_model(test_data, field = FALSE, iterations = 2)

  print_output <- print(result)

  expect_is(print_output, 'disag_model')
  expect_equal(print_output, result)

})

test_that("Check summary.disag_predictions function works as expected", {

  skip_if_not_installed('INLA')
  skip_on_cran()

  result <- disag_model(test_data, iid = FALSE, iterations = 100,
                        list(priormean_intercept = 0,
                             priorsd_intercept = 0.1,
                             priormean_slope = 0.0,
                             priorsd_slope = 0.1,
                             prior_rho_min = 5,
                             prior_rho_prob = 0.01,
                             prior_sigma_max = 0.1,
                             prior_sigma_prob = 0.01,
                             prior_iideffect_sd_max = 0.00001,
                             prior_iideffect_sd_prob = 0.01))

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

  skip_if_not_installed('INLA')
  skip_on_cran()

  result <- disag_model(test_data, iid = FALSE, iterations = 100,
                        list(priormean_intercept = 0,
                             priorsd_intercept = 0.1,
                             priormean_slope = 0.0,
                             priorsd_slope = 0.1,
                             prior_rho_min = 5,
                             prior_rho_prob = 0.01,
                             prior_sigma_max = 0.1,
                             prior_sigma_prob = 0.01,
                             prior_iideffect_sd_max = 0.0001,
                             prior_iideffect_sd_prob = 0.01))

  pred <- predict(result)

  print_output <- print(pred)

  expect_is(print_output, 'disag_prediction')
  expect_equal(print_output, pred)

})
