context("Plotting data")

test_that("Check plot_polygon_data function works as expected", {
  
  polygons <- list()
  for(i in 1:100) {
    row <- ceiling(i/10)
    col <- ifelse(i %% 10 != 0, i %% 10, 10)
    xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
    polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
  }
  
  polys <- do.call(raster::spPolygons, polygons)
  response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  p <- plot_polygon_data(spdf)
  expect_error(plot_polygon_data(polys))
  expect_is(p, 'ggplot')
  
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
  
  polygons <- list()
  for(i in 1:100) {
    row <- ceiling(i/10)
    col <- ifelse(i %% 10 != 0, i %% 10, 10)
    xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
    polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
  }
  
  polys <- do.call(raster::spPolygons, polygons)
  response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  my_mesh <- build_mesh(spdf)
  
  p <- plot_inla_mesh(my_mesh)
  expect_error(plot_inla_mesh(spdf))
  expect_is(p, 'NULL')
  
})

test_that("Check plot.disag.data function works as expected", {
  
  polygons <- list()
  for(i in 1:100) {
    row <- ceiling(i/10)
    col <- ifelse(i %% 10 != 0, i %% 10, 10)
    xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
    polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
  }
  
  polys <- do.call(raster::spPolygons, polygons)
  response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  # Create raster stack
  r <- raster::raster(ncol=20, nrow=20)
  r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
  r2 <- raster::raster(ncol=20, nrow=20)
  r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
  cov_rasters <- raster::stack(r, r2)
  
  result <- prepare_data(polygon_shapefile = spdf, 
                         covariate_rasters = cov_rasters)
  
  p <- plot(result)
  
  expect_error(plot(polys))
  expect_is(p, 'list')
  expect_equal(length(p), 2)
  expect_equal(names(p), c('polygon', 'covariates'))
  
})

test_that("Check plot.fit.result function works as expected", {
  
  fit_result <- get(load(paste0(tempdir(), '/test_fit_result.RData')))
  
  p <- plot(fit_result)
  expect_is(p, 'list')
  
})

