context("Predict model")

polygons <- list()
for(i in 1:100) {
  row <- ceiling(i/10)
  col <- ifelse(i %% 10 != 0, i %% 10, 10)
  xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
  polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
}

polys <- do.call(raster::spPolygons, polygons)
response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 1000))
spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)

# Create raster stack
r <- raster::raster(ncol=20, nrow=20)
r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
r2 <- raster::raster(ncol=20, nrow=20)
r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
cov_stack <- raster::stack(r, r2)

test_data <- prepare_data(polygon_shapefile = spdf, 
                          covariate_rasters = cov_stack)

result <- fit_model(test_data, its = 2)
result_nofield <- fit_model(test_data, its = 2, field = FALSE)


test_that("Check predict_model function works as expected", {
  
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
  
})

test_that("Check predict_uncertainty function works as expected", {
  
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



test_that("Check predict_model function works with newdata", {
  
  newdata <- raster::crop(raster::stack(r, r2), c(0, 180, -90, 90))
  preds1 <- predict_model(result)
  preds2 <- predict_model(result, newdata)
  
  expect_is(preds2, 'predictions')
  expect_equal(length(preds2), 4)
  expect_equal(names(preds2), c('prediction', 'field', 'iid', 'covariates'))
  expect_is(preds2$prediction, 'Raster')
  expect_is(preds2$field, 'Raster')
  expect_true(is.null(preds2$iid))
  expect_is(preds2$covariates, 'Raster')

  expect_false(identical(raster::extent(preds1$prediction), raster::extent(preds2$prediction)))
  
})

test_that("Check predict_uncertainty function works with newdata as expected", {
  
  newdata <- raster::crop(raster::stack(r, r2), c(0, 180, -90, 90))
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


test_that('Check that check_newdata works', {
  
  newdata <- raster::crop(raster::stack(r, r2), c(0, 180, -90, 90))
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
