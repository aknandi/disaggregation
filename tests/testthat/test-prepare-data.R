
context("Preparing data")

test_that("Check prepare_data function works as expected", {
  
  spdf <- rgdal::readOGR(paste0(tempdir(), '/test_spdf.shp'))

  cov_rasters <- raster::stack(paste0(tempdir(), '/test_cov_stack.tif'))
  
  result <- prepare_data(polygon_shapefile = spdf, 
                         covariate_rasters = cov_rasters, 
                         mesh.args = list(max.edge = c(4.0, 8.0)))
  
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
  
  polygon_data <- read.csv(paste0(tempdir(), '/test_polygon_data.csv'))
  cov_data <- read.csv(paste0(tempdir(), '/test_cov_data.csv'))
  coords <- get(load(paste0(tempdir(), '/test_coords.RData')))
  startendindex <- get(load(paste0(tempdir(), '/test_startendindex.RData')))
  mesh <- get(load(paste0(tempdir(), '/test_mesh.RData')))
  
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

