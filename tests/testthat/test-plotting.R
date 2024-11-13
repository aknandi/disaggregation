context("Plotting data")

test_that("Check plot_polygon_data function works as expected", {

  p <- plot_polygon_data(spdf, list(id_var = 'area_id', response_var = 'response'))
  expect_error(plot_polygon_data(polys, list(id_var = 'area_id', response_var = 'response')))
  expect_is(p, 'ggplot')

  p2 <- plot_polygon_data(spdf2, list(id_var = 'area_id', response_var = 'n_positive'))
  expect_is(p2, 'ggplot')

})

test_that("Check plot.disag.data function works as expected", {

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

  fit_result_nofield <- disag_model(test_data, iterations = 100, field = FALSE, family = "poisson", link = "log")

  p1 <- plot(result)

  p2 <- plot(fit_result_nofield)

  p3 <- plot(result, include_iid = TRUE)

  expect_is(p1, 'list')
  expect_equal(length(p1), 2)

  expect_is(p2, 'list')
  expect_equal(length(p2), 2)

  expect_is(p3, 'list')
  expect_equal(length(p3), 2)

})

test_that("Check plot.disag_prediction function works as expected", {

  pred <- predict(result)
  p <- plot(pred)

  expect_is(p, 'gg')

})

