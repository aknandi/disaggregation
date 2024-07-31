
context("Fitting model")

test_that("disag_model produces errors when expected", {

  expect_error(disag_model(list()))
  expect_error(disag_model(test_data, iterations = 'iterations'))
  expect_error(disag_model(test_data, iid = FALSE, priors = list(polygon_sd_men = 0.3, polygon_sd_sd = 0.4)))
  expect_error(disag_model(test_data, priors = c(polygon_sd_mean = 1.2)))
  expect_error(disag_model(test_data, family = 'banana'))
  expect_error(disag_model(test_data, link = 'apple'))

})

test_that("disag_model behaves as expected", {

  result <- disag_model(test_data, iterations = 100, family = 'poisson', link = 'log')

  expect_is(result, 'disag_model')
  expect_equal(length(result), 5)
  expect_equal(length(result$sd_out$par.fixed), terra::nlyr(test_data$covariate_rasters) + 4)
  expect_equal(unique(names(result$sd_out$par.random)), c("iideffect", "nodemean"))
  expect_true(all(c("layer1", "layer2") %in% names(result$sd_out$par.fixed)))
  expect_false(any(names(result$sd_out$par.fixed) == "slope"))
  expect_true(all(c("layer1", "layer2") %in% names(result$opt$par)))
  expect_false(any(names(result$opt$par) == "slope"))

})

test_that("disag_model with 1 covariate behaves as expected", {

  test_data2 <- test_data
  test_data2$covariate_rasters <- test_data2$covariate_rasters[[1]]
  test_data2$covariate_data <- test_data2$covariate_data[, 1:3]

  result <- disag_model(test_data2, iterations = 100, iid = FALSE, family = 'poisson', link = 'log')

  expect_is(result, 'disag_model')
  expect_equal(length(result), 5)

  # Should be intercept, 1 slope, tau gaussian, and 2 for space. None for iid anymore.
  expect_equal(length(result$sd_out$par.fixed), terra::nlyr(test_data2$covariate_rasters) + 3)
  expect_equal(unique(names(result$sd_out$par.random)), c("nodemean"))

  # Confirm only one covariate was fitted.
  expect_equal(sum(names(result$opt$par) == "layer1"), 1)
  expect_false(any(names(result$opt$par) == "layer2"))

})
test_that("user defined model setup is working as expected", {

  binom_data <- prepare_data(polygon_shapefile = spdf_binom,
                             covariate_rasters = cov_stack,
                             sample_size_var = 'sample_size')

  result2 <- disag_model(test_data, iterations = 100, field = FALSE, family = 'poisson', link = 'log')
  result3 <- disag_model(binom_data, iterations = 100, iid = FALSE, family = 'binomial', link = 'logit')
  result4 <- disag_model(test_data, iterations = 100, field = FALSE, iid = FALSE, link = 'identity')

  expect_error(disag_model(test_data, iterations = 100, iid = FALSE, family = 'binomial', link = 'logit'))

  expect_is(result2, 'disag_model')
  expect_equal(length(result2), 5)
  expect_equal(length(result2$sd_out$par.fixed), terra::nlyr(test_data$covariate_rasters) + 2)
  expect_equal(unique(names(result2$sd_out$par.random)), c("iideffect"))
  expect_false(result2$model_setup$field)
  expect_true(result2$model_setup$iid)
  expect_equal(result2$model_setup$family, 'poisson')
  expect_equal(result2$model_setup$link, 'log')

  expect_is(result3, 'disag_model')
  expect_equal(length(result3), 5)
  expect_equal(length(result3$sd_out$par.fixed), terra::nlyr(binom_data$covariate_rasters) + 3)
  expect_equal(unique(names(result3$sd_out$par.random)), c("nodemean"))
  expect_true(result3$model_setup$field)
  expect_false(result3$model_setup$iid)
  expect_equal(result3$model_setup$family, 'binomial')
  expect_equal(result3$model_setup$link, 'logit')

  expect_is(result4, 'disag_model')
  expect_equal(length(result4), 5)
  expect_equal(length(result4$sd_out$par.fixed), terra::nlyr(test_data$covariate_rasters) + 2)
  expect_equal(unique(names(result4$sd_out$par.random)), NULL)
  expect_false(result4$model_setup$field)
  expect_false(result4$model_setup$iid)
  expect_equal(result4$model_setup$family, 'gaussian')
  expect_equal(result4$model_setup$link, 'identity')
})

test_that("make_model_object behaves as expected", {

  result <- make_model_object(test_data, family = 'poisson', link = 'log')

  expect_is(result, 'list')
  expect_equal(sum(sapply(c("par", "fn", "gr", "report"), function(x) !(x %in% names(result)))), 0)

})

test_that("setup_hess_control behaves as expected", {

  obj <- make_model_object(test_data, family = 'poisson', link = 'log')

  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 2, trace = 0))

  hess_control <- setup_hess_control(opt, hess_control_parscale = rep(c(0.9, 1.1), 3), hess_control_ndeps = 1e-3)

  expect_is(hess_control, 'list')
  expect_equal(length(hess_control$parscale), length(opt$par))
  expect_equal(length(hess_control$ndeps), length(opt$par))

})

