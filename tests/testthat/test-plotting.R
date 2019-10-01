context("Plotting data")

polygons <- list()
for(i in 1:100) {
  row <- ceiling(i/10)
  col <- ifelse(i %% 10 != 0, i %% 10, 10)
  xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
  polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
}

polys <- do.call(raster::spPolygons, polygons)
response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 1), sample_size = floor(runif(100, min = 1, max = 100)))
response_df2 <- data.frame(area_id = 1:100, n_positive = runif(100, min = 0, max = 1), sample_size = floor(runif(100, min = 1, max = 100)))
spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
spdf2 <- sp::SpatialPolygonsDataFrame(polys, response_df2)

# Create raster stack
r <- raster::raster(ncol=20, nrow=20)
r <- raster::setExtent(r, raster::extent(spdf))
r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
r2 <- raster::raster(ncol=20, nrow=20)
r2 <- raster::setExtent(r2, raster::extent(spdf))
r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
cov_rasters <- raster::stack(r, r2)

test_data <- prepare_data(polygon_shapefile = spdf, 
                          covariate_rasters = cov_rasters)

test_data2 <- prepare_data(polygon_shapefile = spdf2, 
                           covariate_rasters = cov_rasters,
                           response_var = 'n_positive')

fit_result <- fit_model(test_data, its = 2)

fit_result_nofield <- fit_model(test_data, its = 2, field = FALSE)


test_that("Check plot_polygon_data function works as expected", {
  
  p <- plot_polygon_data(spdf)
  expect_error(plot_polygon_data(polys))
  expect_is(p, 'ggplot')
  
  p2 <- plot_polygon_data(spdf2, zcol = 'n_positive')
  expect_error(plot_polygon_data(spdf2))
  expect_is(p2, 'ggplot')
  
})

test_that("Check plot_covariate_data function works as expected", {
  
  # Create raster stack
  r <- raster::raster(ncol=20, nrow=20)
  r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
  r2 <- raster::raster(ncol=20, nrow=20)
  r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
  cov_stack <- raster::stack(r, r2)
  
  p <- plot_covariate_data(cov_stack)
  expect_error(plot_covariate_data(r))
  expect_is(p, 'trellis')
  
})

test_that("Check plot_inla_mesh function works as expected", {
  
  my_mesh <- build_mesh(spdf)
  
  p <- plot_inla_mesh(my_mesh)
  expect_error(plot_inla_mesh(spdf))
  expect_is(p, 'NULL')
  
})

test_that("Check plot.disag.data function works as expected", {
  
  p <- plot(test_data)
  
  expect_is(p, 'list')
  expect_equal(length(p), 2)
  expect_equal(names(p), c('polygon', 'covariates'))
  
  p2 <- plot(test_data2, zcol = 'n_positive')
  expect_error(plot(test_data2))
  
  expect_is(p2, 'list')
  expect_equal(length(p2), 2)
  expect_equal(names(p2), c('polygon', 'covariates'))
  
})

test_that("Check plot.fit.result function works as expected", {
  
  p1 <- plot(fit_result)
  
  p2 <- plot(fit_result_nofield)
  
  expect_is(p1, 'list')
  expect_equal(length(p1), 2)
  
  expect_is(p2, 'list')
  expect_equal(length(p2), 2)
  
  
})

test_that("Check plot.predictions function works as expected", {
  
  preds <- predict_model(fit_result)
  p1 <- plot(preds)
  
  preds_nofield <- predict_model(fit_result_nofield)
  p2 <- plot(preds_nofield)
  
  preds_withiid <- predict_model(fit_result, predict_iid = TRUE)
  p3 <- plot(preds_withiid)
  
  unc <- predict_uncertainty(fit_result)
  p4 <- plot(unc)
  
  expect_is(p1, 'list')
  expect_equal(length(p1), 3)
  expect_is(p1[[1]], 'trellis')
  expect_is(p1[[2]], 'trellis')
  expect_is(p1[[3]], 'trellis')
  
  expect_is(p2, 'list')
  expect_equal(length(p2), 2)
  expect_is(p2[[1]], 'trellis')
  expect_is(p2[[2]], 'trellis')
  
  expect_is(p3, 'list')
  expect_equal(length(p3), 4)
  expect_is(p3[[1]], 'trellis')
  expect_is(p3[[2]], 'trellis')
  expect_is(p3[[3]], 'trellis')
  expect_is(p3[[4]], 'trellis')
  
  expect_is(p4, 'trellis')
  
})

