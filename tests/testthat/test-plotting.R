context("Plotting data")

test_that("Check plot_polygon_data function works as expected", {

  skip_on_cran()

  p <- plot_polygon_data(spdf, list(id_var = 'area_id', response_var = 'response'))
  expect_error(plot_polygon_data(polys, list(id_var = 'area_id', response_var = 'response')))
  expect_is(p, 'ggplot')

  p2 <- plot_polygon_data(spdf2, list(id_var = 'area_id', response_var = 'n_positive'))
  expect_is(p2, 'ggplot')

})

test_that("Check plot.disag.data function works as expected", {

  skip_if_not_installed('INLA')
  skip_on_cran()

  test_data2 <- prepare_data(polygon_shapefile = spdf2,
                             covariate_rasters = cov_stack,
                             response_var = 'n_positive')

  p <- plot(test_data)

  expect_is(p, 'list')
  expect_equal(length(p), 3)
  expect_equal(names(p), c('polygon', 'covariates', 'mesh'))

  p2 <- plot(test_data2)

  expect_is(p2, 'list')
  expect_equal(length(p2), 3)
  expect_equal(names(p2), c('polygon', 'covariates', 'mesh'))

  p3 <- plot(test_data, which = c(1,3))

  expect_is(p3, 'list')
  expect_equal(length(p3), 2)
  expect_equal(names(p3), c('polygon', 'mesh'))

})

test_that("Check plot.disag_model function works as expected", {

  skip_if_not_installed('INLA')
  skip_on_cran()

  fit_result <- disag_model(test_data, iterations = 10)

  fit_result_nofield <- disag_model(test_data, iterations = 10, field = FALSE)

  p1 <- plot(fit_result)

  p2 <- plot(fit_result_nofield)

  expect_is(p1, 'list')
  expect_equal(length(p1), 2)

  expect_is(p2, 'list')
  expect_equal(length(p2), 2)


})

test_that("Check plot.disag_prediction function works as expected", {

  skip_if_not_installed('INLA')
  skip_on_cran()

  fit_result <- disag_model(test_data, iterations = 1000,
                        iid = TRUE,
                        field = TRUE,
                        family = 'poisson',
                        link = 'log',
                        priors = list(priormean_intercept = 0,
                                      priorsd_intercept = 0.1,
                                      priormean_slope = 0.0,
                                      priorsd_slope = 0.1,
                                      prior_rho_min = 5,
                                      prior_rho_prob = 0.01,
                                      prior_sigma_max = 0.1,
                                      prior_sigma_prob = 0.01,
                                      prior_iideffect_sd_max = 0.0001,
                                      prior_iideffect_sd_prob = 0.01))
  pred <- predict(fit_result)
  p <- plot(pred)

  expect_is(p, 'gg')

})

