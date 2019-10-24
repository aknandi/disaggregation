context("Predict model")

polygons <- list()
n_polygon_per_side <- 10
n_polygons <- n_polygon_per_side * n_polygon_per_side
n_pixels_per_side <- n_polygon_per_side * 2

for(i in 1:n_polygons) {
  row <- ceiling(i/n_polygon_per_side)
  col <- ifelse(i %% n_polygon_per_side != 0, i %% n_polygon_per_side, n_polygon_per_side)
  xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
  polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
}

polys <- do.call(raster::spPolygons, polygons)
response_df <- data.frame(area_id = 1:n_polygons, response = runif(n_polygons, min = 0, max = 1000))
spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)

# Create raster stack
r <- raster::raster(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
r <- raster::setExtent(r, raster::extent(spdf))
r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side), 3))
r2 <- raster::raster(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
r2 <- raster::setExtent(r2, raster::extent(spdf))
r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/n_pixels_per_side), 3))
cov_stack <- raster::stack(r, r2)

test_data <- prepare_data(polygon_shapefile = spdf, 
                          covariate_rasters = cov_stack)


test_that("Check predict_model and predict_uncertainty function works as expected", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  result <- fit_model(test_data, iterations = 2)
  
  result_nofield <- fit_model(test_data, iterations = 2, field = FALSE)
  
  preds <- predict_model(result)
  
  expect_is(preds, 'predictions')
  expect_equal(length(preds), 4)
  expect_equal(names(preds), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(preds$prediction, 'Raster')
  expect_is(preds$field, 'Raster')
  expect_true(is.null(preds$iid))
  expect_is(preds$covariates, 'Raster')
  
  
  preds2 <- predict_model(result_nofield)
  
  expect_is(preds2, 'predictions')
  expect_equal(length(preds2), 4)
  expect_equal(names(preds2), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(preds2$prediction, 'Raster')
  expect_true(is.null(preds2$field))
  expect_true(is.null(preds2$iid))
  expect_is(preds2$covariates, 'Raster')
  
  
  preds3 <- predict_model(result, predict_iid = TRUE)
  
  expect_is(preds3, 'predictions')
  expect_equal(length(preds3), 4)
  expect_equal(names(preds3), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(preds3$prediction, 'Raster')
  expect_is(preds3$field, 'Raster')
  expect_is(preds3$iid, 'Raster')
  expect_is(preds3$covariates, 'Raster')
  

  unc <- predict_uncertainty(result)
  
  expect_is(unc, 'uncertainty')
  expect_equal(length(unc), 2)
  expect_equal(names(unc), c('realisations', 'predictions_ci'))
  expect_is(unc$realisations, 'RasterStack')
  expect_is(unc$predictions_ci, 'RasterBrick')
  expect_equal(raster::nlayers(unc$realisations), 100)
  expect_equal(raster::nlayers(unc$predictions_ci), 2)
  
  
  unc2 <- predict_uncertainty(result_nofield, N = 10)
  
  expect_is(unc2, 'uncertainty')
  expect_equal(length(unc2), 2)
  expect_equal(names(unc2), c('realisations', 'predictions_ci'))
  expect_is(unc2$realisations, 'RasterStack')
  expect_is(unc2$predictions_ci, 'RasterBrick')
  expect_equal(raster::nlayers(unc2$realisations), 10)
  expect_equal(raster::nlayers(unc2$predictions_ci), 2)
  
  
  unc2 <- predict_uncertainty(result, predict_iid = TRUE, N = 10)
  
  expect_is(unc2, 'uncertainty')
  expect_equal(length(unc2), 2)
  expect_equal(names(unc2), c('realisations', 'predictions_ci'))
  expect_is(unc2$realisations, 'RasterStack')
  expect_is(unc2$predictions_ci, 'RasterBrick')
  expect_equal(raster::nlayers(unc2$realisations), 10)
  expect_equal(raster::nlayers(unc2$predictions_ci), 2)
  
})



test_that("Check predict_model and predict_uncertainty function works with newdata", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  result <- fit_model(test_data, field = FALSE, iid = TRUE, iterations = 2)
  
  newdata <- raster::crop(raster::stack(r, r2), c(0, 10, 0, 10))
  preds1 <- predict_model(result)
  preds2 <- predict_model(result, newdata, predict_iid = TRUE)
  
  expect_is(preds2, 'predictions')
  expect_equal(length(preds2), 4)
  expect_equal(names(preds2), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(preds2$prediction, 'Raster')
  expect_true(is.null(preds2$field))
  expect_is(preds2$iid, 'Raster')
  expect_is(preds2$covariates, 'Raster')

  expect_false(identical(raster::extent(preds1$prediction), raster::extent(preds2$prediction)))
  
  
  unc1 <- predict_uncertainty(result, N = 5)
  unc2 <- predict_uncertainty(result, newdata, N = 5)
  
  expect_is(unc1, 'uncertainty')
  expect_equal(length(unc2), 2)
  expect_equal(names(unc2), c('realisations', 'predictions_ci'))
  expect_is(unc2$realisations, 'RasterStack')
  expect_is(unc2$predictions_ci, 'RasterBrick')
  expect_equal(raster::nlayers(unc2$realisations), 5)
  expect_equal(raster::nlayers(unc2$predictions_ci), 2)

  expect_false(identical(raster::extent(unc1$realisations), raster::extent(unc2$realisations)))
  
})

test_that('Check that predict.fit.model works', {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  result <- fit_model(test_data, field = FALSE, iterations = 2)
  
  preds <- predict(result, N = 5)
  
  expect_is(preds, 'list')
  expect_equal(length(preds), 2)
  expect_is(preds$mean_predictions, 'predictions')
  expect_is(preds$uncertainty_predictions, 'uncertainty')
  expect_equal(length(preds$mean_predictions), 4)
  expect_equal(names(preds$mean_predictions), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(preds$mean_predictions$prediction, 'Raster')
  expect_true(is.null(preds$mean_predictions$field))
  expect_true(is.null(preds$mean_predictions$iid))
  expect_is(preds$mean_predictions$covariates, 'Raster')
  expect_equal(length(preds$uncertainty_predictions), 2)
  expect_equal(names(preds$uncertainty_predictions), c('realisations', 'predictions_ci'))
  expect_equal(raster::nlayers(preds$uncertainty_predictions$realisations), 5)
  expect_equal(raster::nlayers(preds$uncertainty_predictions$predictions_ci), 2)

  
})


test_that('Check that check_newdata works', {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  result <- fit_model(test_data, field = FALSE, iterations = 2)
  
  newdata <- raster::crop(raster::stack(r, r2), c(0, 10, 0, 10))
  nd1 <- check_newdata(newdata, result)
  expect_is(nd1, 'RasterBrick')
  
  nn <- newdata[[1]]
  names(nn) <- 'extra_uneeded'
  newdata2 <- raster::stack(newdata, nn)
  expect_error(check_newdata(newdata2, result), NA)

  newdata3 <- newdata[[1]]
  expect_error(check_newdata(newdata3, result), 'All covariates')
  
  newdata4 <- result$data$covariate_data
  expect_error(check_newdata(newdata4, result), 'newdata should be NULL or')
  
  
})

test_that('Check that setup_objects works', {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  result <- fit_model(test_data, iterations = 2)
  
  objects <- setup_objects(result)
  
  expect_is(objects, 'list')
  expect_equal(length(objects), 3)
  expect_equal(names(objects), c('covariates', 'field_objects', 'iid_objects'))
  expect_is(objects$field_objects, 'list')
  expect_true(is.null(objects$iid_objects))

  newdata <- raster::crop(raster::stack(r, r2), c(0, 180, -90, 90))
  objects2 <- setup_objects(result, newdata)
  
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
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  result <- fit_model(test_data, iterations = 2)
  
  objects <- setup_objects(result)
  
  pars <- result$obj$env$last.par.best
  pars <- split(pars, names(pars))
  
  pred <- predict_single_raster(pars, 
                                objects = objects,
                                link_function = result$model_setup$link)
  
  expect_is(pred, 'list')
  expect_equal(length(pred), 4)
  expect_equal(names(pred), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(pred$prediction, 'Raster')
  expect_is(pred$field, 'Raster')
  expect_true(is.null(pred$iid))
  expect_is(pred$covariates, 'Raster')
  
  objects2 <- setup_objects(result, predict_iid = TRUE)
  
  pred <- predict_single_raster(pars, 
                                objects = objects2,
                                link_function = result$model_setup$link)
  
  expect_is(pred, 'list')
  expect_equal(length(pred), 4)
  expect_equal(names(pred), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(pred$prediction, 'Raster')
  expect_is(pred$field, 'Raster')
  expect_is(pred$iid, 'Raster')
  expect_is(pred$covariates, 'Raster')
  
})
