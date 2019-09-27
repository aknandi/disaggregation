context("Summary functions")

polygons <- list()
for(i in 1:100) {
  row <- ceiling(i/10)
  col <- ifelse(i %% 10 != 0, i %% 10, 10)
  xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
  polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
}

polys <- do.call(raster::spPolygons, polygons)
response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 1), sample_size = floor(runif(100, min = 1, max = 100)))
spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)

# Create raster stack
r <- raster::raster(ncol=20, nrow=20)
r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
r2 <- raster::raster(ncol=20, nrow=20)
r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
cov_stack <- raster::stack(r, r2)

test_data <- prepare_data(polygon_shapefile = spdf, 
                          covariate_rasters = cov_stack)


test_that("Check summary.disag.data function works as expected", {
  
  data_summary <- summary(test_data)
  
  expect_is(data_summary, 'list')
  expect_equal(length(data_summary), 3)
  expect_equal(names(data_summary), c('number_polygons', 'number_covariates', 'covariate_summary'))
  expect_is(data_summary$number_polygons, 'integer')
  expect_is(data_summary$number_covariates, 'integer')
  expect_is(data_summary$covariate_summary, 'table')
  expect_equal(ncol(data_summary$covariate_summary), data_summary$number_covariates)
  
})

test_that("Check summary.fit.model function works as expected", {
  
  result <- fit_model(test_data, its = 2)
  
  model_summary <- summary(result)
  
  expect_is(model_summary, 'list')
  expect_equal(length(model_summary), 3)
  expect_equal(names(model_summary), c('model_params', 'nll', 'metrics'))
  expect_is(model_summary$model_params, 'matrix')
  expect_is(model_summary$nll, 'numeric')
  expect_is(model_summary$metrics, 'data.frame')
  expect_equal(dim(model_summary$metrics), c(1, 5))
  
})