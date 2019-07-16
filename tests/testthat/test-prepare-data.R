
context("Preparing data")

test_that("Check prepare_data function works as expected", {
  
  cds1 <- rbind(c(-120,-20), c(-100,5), c(-60, 0), c(-100,-60), c(-120,-20))
  cds2 <- rbind(c(30,0), c(50,60), c(70,0), c(70,-55), c(30,0))
  cds3 <- rbind(c(10,20), c(10,50), c(-60,40), c(-30,-30), c(30,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('1', '2', '3'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)

  # Create raster stack
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  cov_rasters <- raster::stack(r, r)
  
  result <- prepare_data(polygon_shapefile = spdf, 
                         covariate_rasters = cov_rasters)
  
  expect_is(result, 'disag.data')
  expect_equal(length(result), 5)
  expect_equal(names(result), c('polygon_data', 'covariate_data', 'coords', 'startendindex', 'mesh'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$coords, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coords))
  
})

test_that("Check as.disag.data function works as expected", {
  
  cds1 <- rbind(c(-120,-20), c(-100,5), c(-60, 0), c(-100,-60), c(-120,-20))
  cds2 <- rbind(c(30,0), c(50,60), c(70,0), c(70,-55), c(30,0))
  cds3 <- rbind(c(10,20), c(10,50), c(-60,40), c(-30,-30), c(30,10))
  polys <- raster::spPolygons(cds1, cds2, cds3)
  
  response_df <- data.frame(area_id = c('1', '2', '3'), response = c(4, 7, 2))
  
  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
  
  polygon_data <- getPolygonData(spdf, id_var = 'area_id', response_var = 'response')
  
  # Create raster stack
  r <- raster::raster(ncol=36, nrow=18)
  r[] <- 1:raster::ncell(r)
  r <- raster::stack(r, r)
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  cov_data <- parallelExtract(r, spdf, fun = NULL, id = 'area_id')
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  coords <- extractCoordsForMesh(r, cov_data)
  
  startendindex <- getStartendindex(cov_data, polygon_data, 'area_id')
  
  mesh <- build_mesh(spdf)

  result <- as.disag.data(polygon_data, cov_data, coords, startendindex, mesh)
  
  expect_is(result, 'disag.data')
  expect_equal(length(result), 5)
  expect_equal(names(result), c('polygon_data', 'covariate_data', 'coords', 'startendindex', 'mesh'))
  expect_is(result$polygon_data, 'data.frame')
  expect_is(result$covariate_data, 'data.frame')
  expect_is(result$coords, 'matrix')
  expect_is(result$startendindex, 'matrix')
  expect_is(result$mesh, 'inla.mesh')
  expect_equal(nrow(result$polygon_data), nrow(result$startendindex))
  expect_equal(nrow(result$covariate_data), nrow(result$coords))
  
})

