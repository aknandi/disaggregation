
context("Preparing data")

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
cov_rasters <- raster::stack(r, r2)


test_that("Check prepare_data function works as expected", {
  
  result <- prepare_data(polygon_shapefile = spdf, 
                         covariate_rasters = cov_rasters)
  
  expect_is(result, 'disag.data')
  expect_equal(length(result), 8)
  expect_equal(names(result), c('polygon_shapefile', 'covariate_rasters', 'polygon_data', 'covariate_data', 'aggregation_pixels', 'coords', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coords, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(sum(is.na(result$polygon_data$N)), length(result$polygon_data$N))
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coords))
  
})

test_that("Check prepare_data function with sample size works as expected", {
  
  result <- prepare_data(polygon_shapefile = spdf, 
                         covariate_rasters = cov_rasters,
                         sample_size_var = 'sample_size')
  
  expect_is(result, 'disag.data')
  expect_equal(length(result), 8)
  expect_equal(names(result), c('polygon_shapefile', 'covariate_rasters', 'polygon_data', 'covariate_data', 'aggregation_pixels', 'coords', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coords, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(sum(is.na(result$polygon_data$N)), 0)
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coords))
  
})


test_that("Check as.disag.data function works as expected", {
  

  polygon_data <- getPolygonData(spdf, id_var = 'area_id', response_var = 'response')
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  cov_data <- parallelExtract(cov_rasters, spdf, fun = NULL, id = 'area_id')
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  aggregation_data <- rep(1, nrow(cov_data))
  
  coords <- extractCoordsForMesh(cov_rasters, cov_data)
  
  startendindex <- getStartendindex(cov_data, polygon_data, 'area_id')
  
  mesh <- build_mesh(spdf)

  result <- as.disag.data(spdf, 
                          cov_rasters,
                          polygon_data, 
                          cov_data, 
                          aggregation_data,
                          coords, 
                          startendindex, 
                          mesh)
  
  expect_is(result, 'disag.data')
  expect_equal(length(result), 8)
  expect_equal(names(result), c('polygon_shapefile', 'covariate_rasters', 'polygon_data', 'covariate_data', 'aggregation_pixels', 'coords', 'startendindex', 'mesh'))
  expect_is(result$polygon_shapefile, 'SpatialPolygonsDataFrame')
  expect_is(result$covariate_rasters, c('RasterBrick', 'RasterStack'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$aggregation_pixels, 'numeric')
  expect_is(result$coords, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coords))
  
})

