
context("Fitting model")

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
response_binom_df <- data.frame(area_id = 1:n_polygons, response = N*runif(n_polygons, min = 0, max = 1), sample_size = N)

spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
spdf_binom <- sp::SpatialPolygonsDataFrame(polys, response_binom_df)

# Create raster stack
r <- raster::raster(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
r <- raster::setExtent(r, raster::extent(spdf))
r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side), 3))
r2 <- raster::raster(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
r2 <- raster::setExtent(r2, raster::extent(spdf))
r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/n_pixels_per_side), 3))
cov_stack <- raster::stack(r, r2)

if(identical(Sys.getenv("NOT_CRAN"), "true")) {
  test_data <- prepare_data(polygon_shapefile = spdf, 
                            covariate_rasters = cov_stack)
} else {
  test_data <- prepare_data(polygon_shapefile = spdf, 
                            covariate_rasters = cov_stack,
                            makeMesh = FALSE)
}

test_that("disag_model produces errors when expected", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  expect_error(disag_model(list()))
  expect_error(disag_model(test_data, iterations = 'iterations'))
  expect_error(disag_model(test_data, priors = list(polygon_sd_men = 0.3, polygon_sd_sd = 0.4)))
  expect_error(disag_model(test_data, priors = c(polygon_sd_mean = 1.2)))
  expect_error(disag_model(test_data, family = 'banana'))
  expect_error(disag_model(test_data, link = 'apple'))
  
})

test_that("disag_model behaves as expected", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  result <- disag_model(test_data, iterations = 2, iid = FALSE)

  expect_is(result, 'disag_model')
  expect_equal(length(result), 5)
  expect_equal(length(result$sd_out$par.fixed), raster::nlayers(test_data$covariate_rasters) + 4)
  expect_equal(unique(names(result$sd_out$par.random)), c("nodemean"))
  
  
  
})




test_that("disag_model with 1 covariate behaves as expected", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  test_data2 <- test_data
  test_data2$covariate_rasters <- test_data2$covariate_rasters[[1]]
  test_data2$covariate_data <- test_data2$covariate_data[, 1:3]
  
  result <- disag_model(test_data2, iterations = 2, iid = FALSE)
  
  expect_is(result, 'disag_model')
  expect_equal(length(result), 5)
  
  # Should be intercept, 1 slope, tau gaussian, and 2 for space. None for iid anymore.
  expect_equal(length(result$sd_out$par.fixed), raster::nlayers(test_data$covariate_rasters) + 3)
  expect_equal(unique(names(result$sd_out$par.random)), c("nodemean"))
  
  # Confirm only two covariates were fitted.
  expect_equal(sum(names(result$opt$par) == 'slope'), 1)
  
})
test_that("user defined model setup is working as expected", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  binom_data <- prepare_data(polygon_shapefile = spdf_binom, 
                             covariate_rasters = cov_stack,
                             sample_size_var = 'sample_size')
  
  result2 <- disag_model(test_data, iterations = 2, field = FALSE, family = 'poisson', link = 'log')
  result3 <- disag_model(binom_data, iterations = 2, iid = FALSE, family = 'binomial', link = 'logit')
  result4 <- disag_model(test_data, iterations = 2, field = FALSE, iid = FALSE, link = 'identity')
  
  expect_error(disag_model(test_data, iterations = 2, iid = FALSE, family = 'binomial', link = 'logit'))
  
  expect_is(result2, 'disag_model')
  expect_equal(length(result2), 5)
  expect_equal(length(result2$sd_out$par.fixed), raster::nlayers(test_data$covariate_rasters) + 2)
  expect_equal(unique(names(result2$sd_out$par.random)), c("iideffect"))
  expect_false(result2$model_setup$field)
  expect_true(result2$model_setup$iid)
  expect_equal(result2$model_setup$family, 'poisson')
  expect_equal(result2$model_setup$link, 'log')
  
  expect_is(result3, 'disag_model')
  expect_equal(length(result3), 5)
  expect_equal(length(result3$sd_out$par.fixed), raster::nlayers(binom_data$covariate_rasters) + 3)
  expect_equal(unique(names(result3$sd_out$par.random)), c("nodemean"))
  expect_true(result3$model_setup$field)
  expect_false(result3$model_setup$iid)
  expect_equal(result3$model_setup$family, 'binomial')
  expect_equal(result3$model_setup$link, 'logit')
  
  expect_is(result4, 'disag_model')
  expect_equal(length(result4), 5)
  expect_equal(length(result4$sd_out$par.fixed), raster::nlayers(test_data$covariate_rasters) + 2)
  expect_equal(unique(names(result4$sd_out$par.random)), NULL)
  expect_false(result4$model_setup$field)
  expect_false(result4$model_setup$iid)
  expect_equal(result4$model_setup$family, 'gaussian')
  expect_equal(result4$model_setup$link, 'identity')
})

test_that("make_model_object behaves as expected", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  result <- make_model_object(test_data)
  
  expect_is(result, 'list')
  expect_equal(sum(sapply(c("par", "fn", "gr", "report"), function(x) !(x %in% names(result)))), 0)
  
})

test_that("setup_hess_control behaves as expected", {
  
  skip_if_not_installed('INLA')
  skip_on_cran()
  
  obj <- make_model_object(test_data)
  
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 2, trace = 0))
  
  hess_control <- setup_hess_control(opt, hess_control_parscale = c(rep(c(0.9, 1.1), 3), 1), hess_control_ndeps = 1e-3)
  
  expect_is(hess_control, 'list')
  expect_equal(length(hess_control$parscale), length(opt$par))
  expect_equal(length(hess_control$ndeps), length(opt$par))
  
})

