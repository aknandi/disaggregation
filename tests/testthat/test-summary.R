context("Summary functions")

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

test_that("Check summary.disag.data function works as expected", {
  
  skip_on_cran()
  
  test_data <- prepare_data(polygon_shapefile = spdf, 
                            covariate_rasters = cov_stack)
  
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
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  test_data <- prepare_data(polygon_shapefile = spdf, 
                            covariate_rasters = cov_stack)
  
  result <- fit_model(test_data, field = FALSE, iterations = 2)
  
  model_summary <- summary(result)
  
  expect_is(model_summary, 'list')
  expect_equal(length(model_summary), 3)
  expect_equal(names(model_summary), c('model_params', 'nll', 'metrics'))
  expect_is(model_summary$model_params, 'matrix')
  expect_is(model_summary$nll, 'numeric')
  expect_is(model_summary$metrics, 'data.frame')
  expect_equal(dim(model_summary$metrics), c(1, 5))
  
})