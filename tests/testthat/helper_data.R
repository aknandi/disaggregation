polygons <- list()
n_polygon_per_side <- 10
n_polygons <- n_polygon_per_side * n_polygon_per_side
n_pixels_per_side <- n_polygon_per_side * 2

for(i in 1:n_polygons) {
  row <- ceiling(i/n_polygon_per_side)
  col <- ifelse(i %% n_polygon_per_side != 0, i %% n_polygon_per_side, n_polygon_per_side)
  xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
  polygons[[i]] <- list(cbind(c(xmin, xmax, xmax, xmin, xmin),
                              c(ymax, ymax, ymin, ymin, ymax)))
}

polys <- lapply(polygons, sf::st_polygon)
N <- floor(runif(n_polygons, min = 1, max = 100))
response_df <- data.frame(area_id = 1:n_polygons, response = runif(n_polygons, min = 0, max = 1000))
response_na_df <- data.frame(area_id = 1:n_polygons, response = c(runif(n_polygons - 1, min = 0, max = 1000), NA))
response_binom_df <- data.frame(area_id = 1:n_polygons, response = N*runif(n_polygons, min = 0, max = 1), sample_size = N)
response_df2 <- data.frame(area_id = 1:n_polygons, n_positive = runif(n_polygons, min = 0, max = 1), sample_size = floor(runif(n_polygons, min = 1, max = 100)))

spdf <- sf::st_sf(response_df, geometry = polys)
spdf_na <- sf::st_sf(response_na_df, geometry = polys)
spdf_binom <- sf::st_sf(response_binom_df, geometry = polys)
spdf2 <- sf::st_sf(response_df2, geometry = polys)

# Create raster stack
r <- terra::rast(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
terra::ext(r) <- terra::ext(spdf)
r[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side), 3))
r2 <- terra::rast(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
terra::ext(r2) <- terra::ext(spdf)
r2[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ceiling(x/n_pixels_per_side), 3))
cov_stack <- c(r, r2)
names(cov_stack) <- c('layer1', 'layer2')

test_data <- prepare_data(polygon_shapefile = spdf,
                          covariate_rasters = cov_stack)
