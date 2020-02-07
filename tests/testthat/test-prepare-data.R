
context("Preparing data")

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
N <- floor(runif(n_polygons, min = 1, max = 100))
response_df <- data.frame(area_id = 1:n_polygons, response = runif(n_polygons, min = 0, max = 1000))
response_na_df <- data.frame(area_id = 1:n_polygons, response = c(runif(n_polygons - 1, min = 0, max = 1000), NA))
response_binom_df <- data.frame(area_id = 1:n_polygons, response = N*runif(n_polygons, min = 0, max = 1), sample_size = N)

spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
spdf_na <- sp::SpatialPolygonsDataFrame(polys, response_na_df)
spdf_binom <- sp::SpatialPolygonsDataFrame(polys, response_binom_df)

# Create raster stack
r <- raster::raster(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
r <- raster::setExtent(r, raster::extent(spdf))
r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side), 3))
r2 <- raster::raster(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
r2 <- raster::setExtent(r2, raster::extent(spdf))
r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/n_pixels_per_side), 3))
cov_stack <- raster::stack(r, r2)


test_that("Check prepare_data function works as expected", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  result <- prepare_data(polygon_shapefile = spdf, 
                         covariate_rasters = cov_stack)
  
  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data', 
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$coordsForPrediction, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(sum(is.na(result$polygon_data$N)), length(result$polygon_data$N))
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))
  
})

test_that("Check prepare_data function with sample size works as expected", {
  
  skip_on_cran()
  
  result <- prepare_data(polygon_shapefile = spdf_binom, 
                         covariate_rasters = cov_stack,
                         sample_size_var = 'sample_size',
                         makeMesh = FALSE)
  
  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data', 
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$coordsForPrediction, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_true(is.null(result$mesh))
  expect_equal(sum(is.na(result$polygon_data$N)), 0)
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))
  
})

test_that("Check prepare_data function deals with NAs as expected", {
  
  skip_on_cran()
  
  cov_stack_na <- cov_stack
  cov_stack_na[[1]][c(1:10)] <- NA
  
  aggregation_raster_na <- r
  aggregation_raster_na[c(1:10)] <- NA
  
  expect_error(prepare_data(polygon_shapefile = spdf_na, covariate_rasters = cov_stack, makeMesh = FALSE))
  expect_error(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_stack_na, makeMesh = FALSE))
  expect_error(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_stack, aggregation_raster = aggregation_raster_na, makeMesh = FALSE))
               
  result <- prepare_data(polygon_shapefile = spdf_na, 
                         covariate_rasters = cov_stack_na,
                         aggregation_raster = aggregation_raster_na,
                         na.action = TRUE,
                         makeMesh = FALSE)
  
  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data', 
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_true(is.null(result$mesh))
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))
  expect_equal(sum(is.na(result$polygon_data$response)), 0)
  expect_equal(sum(is.na(result$covariate_data)), 0)
  expect_equal(sum(is.na(result$aggregation_pixels)), 0)
  expect_equal(nrow(result$polygon_shapefile), nrow(spdf_na) - 1)
})


test_that("Check as.disag_data function works as expected", {
  
  skip_on_cran()
  
  polygon_data <- getPolygonData(spdf, id_var = 'area_id', response_var = 'response')
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  cov_data <- parallelExtract(cov_stack, spdf, fun = NULL, id = 'area_id')
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  aggregation_data <- rep(1, nrow(cov_data))
  
  coordsForFit <- extractCoordsForMesh(cov_stack, cov_data$cellid)
  
  coordsForPrediction <- extractCoordsForMesh(cov_stack)
  
  startendindex <- getStartendindex(cov_data, polygon_data, 'area_id')
  
  result <- as.disag_data(spdf, 
                          list('area_id', 'response'),
                          cov_stack,
                          polygon_data, 
                          cov_data, 
                          aggregation_data,
                          coordsForFit, 
                          coordsForPrediction,
                          startendindex,
                          mesh = NULL)
  
  expect_is(result, 'disag_data')
  expect_equal(length(result), 10)
  expect_equal(names(result), c('polygon_shapefile', 'shapefile_names', 'covariate_rasters', 'polygon_data', 'covariate_data', 
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$shapefile_names, 'list')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$coordsForPrediction, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_true(is.null(result$mesh))
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))
  
})

