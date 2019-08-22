
context("Fitting model")

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
cov_stack <- raster::stack(r, r2)

test_data <- prepare_data(polygon_shapefile = spdf, 
                          covariate_rasters = cov_stack)

test_that("fit_model produces errors whe expected", {
  
  expect_error(fit_model(list()))
  expect_error(fit_model(test_data, its = 'its'))
  expect_error(fit_model(test_data, priors = list(polygon_sd_men = 0.3, polygon_sd_sd = 0.4)))
  expect_error(fit_model(test_data, priors = c(polygon_sd_mean = 1.2)))
  expect_is(fit_model(test_data, priors = list(polygon_sd_mean = 0.3, polygon_sd_mean = 0.2)), 'fit.result')
  
})

test_that("fit_model behaves as expected", {
  
  result <- fit_model(test_data, its = 2)

  save(result, file = paste0(tempdir(), '/test_fit_result.RData'))
  
  expect_is(result, 'fit.result')
  expect_equal(length(result), 3)
  expect_equal(length(result$sd_out$par.fixed), raster::nlayers(test_data$covariate_rasters) + 4)
  expect_equal(unique(names(result$sd_out$par.random)), c("iideffect", "nodemean"))
  
})

test_that("field and iid flags in fit_model are working as expected", {
  
  result2 <- fit_model(test_data, its = 2, field = FALSE)
  result3 <- fit_model(test_data, its = 2, iid = FALSE)
  result4 <- fit_model(test_data, its = 2, field = FALSE, iid = FALSE)
  
  expect_is(result2, 'fit.result')
  expect_equal(length(result2), 3)
  expect_equal(length(result2$sd_out$par.fixed), raster::nlayers(test_data$covariate_rasters) + 2)
  expect_equal(unique(names(result2$sd_out$par.random)), c("iideffect"))
  
  expect_is(result3, 'fit.result')
  expect_equal(length(result3), 3)
  expect_equal(length(result3$sd_out$par.fixed), raster::nlayers(test_data$covariate_rasters) + 4)
  expect_equal(unique(names(result3$sd_out$par.random)), c("nodemean"))
  
  expect_is(result4, 'fit.result')
  expect_equal(length(result4), 3)
  expect_equal(length(result4$sd_out$par.fixed), raster::nlayers(test_data$covariate_rasters) + 2)
  expect_equal(unique(names(result4$sd_out$par.random)), NULL)
})