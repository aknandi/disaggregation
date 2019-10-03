
context("Preparing data")

polygons <- list()
for(i in 1:100) {
  row <- ceiling(i/10)
  col <- ifelse(i %% 10 != 0, i %% 10, 10)
  xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
  polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
}

polys <- do.call(raster::spPolygons, polygons)
N <- floor(runif(100, min = 1, max = 100))
response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 1000))
response_na_df <- data.frame(area_id = 1:100, response = c(runif(99, min = 0, max = 1000), NA))
response_binom_df <- data.frame(area_id = 1:100, response = N*runif(100, min = 0, max = 1), sample_size = N)

spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
spdf_na <- sp::SpatialPolygonsDataFrame(polys, response_na_df)
spdf_binom <- sp::SpatialPolygonsDataFrame(polys, response_binom_df)

# Create raster stack
r <- raster::raster(ncol=20, nrow=20)
r <- raster::setExtent(r, raster::extent(spdf))
r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
r2 <- raster::raster(ncol=20, nrow=20)
r2 <- raster::setExtent(r2, raster::extent(spdf))
r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
cov_rasters <- raster::stack(r, r2)


test_that("Check prepare_data function works as expected", {
  
  result <- prepare_data(polygon_shapefile = spdf, 
                         covariate_rasters = cov_rasters)
  
  expect_is(result, 'disag.data')
  expect_equal(length(result), 9)
  expect_equal(names(result), c('polygon_shapefile', 'covariate_rasters', 'polygon_data', 'covariate_data', 
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
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
  
  result <- prepare_data(polygon_shapefile = spdf_binom, 
                         covariate_rasters = cov_rasters,
                         sample_size_var = 'sample_size')
  
  expect_is(result, 'disag.data')
  expect_equal(length(result), 9)
  expect_equal(names(result), c('polygon_shapefile', 'covariate_rasters', 'polygon_data', 'covariate_data', 
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$coordsForPrediction, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(sum(is.na(result$polygon_data$N)), 0)
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))
  
})

test_that("Check prepare_data function deals with NAs as expected", {
  
  cov_rasters_na <- cov_rasters
  cov_rasters_na[[1]][c(1:10)] <- NA
  
  aggregation_raster_na <- r
  aggregation_raster_na[c(1:10)] <- NA
  
  expect_error(prepare_data(polygon_shapefile = spdf_na, covariate_rasters = cov_rasters))
  expect_error(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_rasters_na))
  expect_error(prepare_data(polygon_shapefile = spdf, covariate_rasters = cov_rasters, aggregation_raster = aggregation_raster_na))
               
  result <- prepare_data(polygon_shapefile = spdf_na, 
                         covariate_rasters = cov_rasters_na,
                         aggregation_raster = aggregation_raster_na,
                         na.action = TRUE)
  
  expect_is(result, 'disag.data')
  expect_equal(length(result), 9)
  expect_equal(names(result), c('polygon_shapefile', 'covariate_rasters', 'polygon_data', 'covariate_data', 
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))
  expect_equal(sum(is.na(result$polygon_data$response)), 0)
  expect_equal(sum(is.na(result$covariate_data)), 0)
  expect_equal(sum(is.na(result$aggregation_pixels)), 0)
  expect_equal(nrow(result$polygon_shapefile), nrow(spdf_na) - 1)
})


test_that("Check as.disag.data function works as expected", {
  

  polygon_data <- getPolygonData(spdf, id_var = 'area_id', response_var = 'response')
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  cov_data <- parallelExtract(cov_rasters, spdf, fun = NULL, id = 'area_id')
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  aggregation_data <- rep(1, nrow(cov_data))
  
  coordsForFit <- extractCoordsForMesh(cov_rasters, cov_data$cellid)
  
  coordsForPrediction <- extractCoordsForMesh(cov_rasters)
  
  startendindex <- getStartendindex(cov_data, polygon_data, 'area_id')
  
  mesh <- build_mesh(spdf)

  result <- as.disag.data(spdf, 
                          cov_rasters,
                          polygon_data, 
                          cov_data, 
                          aggregation_data,
                          coordsForFit, 
                          coordsForPrediction,
                          startendindex, 
                          mesh)
  
  expect_is(result, 'disag.data')
  expect_equal(length(result), 9)
  expect_equal(names(result), c('polygon_shapefile', 'covariate_rasters', 'polygon_data', 'covariate_data', 
                                'aggregation_pixels', 'coordsForFit', 'coordsForPrediction', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coordsForFit, 'matrix')
  expect_is(result$coordsForPrediction, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coordsForFit))
  
})

