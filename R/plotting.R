#' These are plot methods for input data and predictions for disaggragation
#'
#' @param x Object to be plotted
#' @param ... Further arguments passed to or from other methods.
#' 
#' @import ggplot2
#' 
#' @name plot.disag.data
#' 
#' @export

plot.disag.data <- function(x, ...) {
  
  # Plot polygon data, covariate rasters and mesh
  plots <- list()
  
  plots$polygon <- plot_polygon_data(x$polygon_shapefile)
  plots$covariates <- plot_covariate_data(x$covariate_rasters)
  plots$mesh <- plot_inla_mesh(x$mesh)
  
  #print(cowplot::plot_grid(plotlist = plots))
  
  return(invisible(plots))
}

#' Plot polygon data from SpatialPolygonDataFrame
#'
#' @param x Object to be plotted
#' @name plot_polygon_data

plot_polygon_data <- function(x) {
  
  stopifnot(inherits(x, 'SpatialPolygonsDataFrame'))
  
  df_fortify <- fortify(x, region = 'area_id')
  
  df <- x@data
  df <- dplyr::mutate(df, area_id = as.character(area_id)) 
  df <- dplyr::left_join(df_fortify, df, by = c('id' = 'area_id'))
  
  p <- ggplot(df, aes(long, lat, group = group, fill = response)) + 
    geom_polygon() +
    coord_equal() +
    scale_fill_viridis_c(trans = 'identity')
  
  print(p)
  return(p)
}

#' Plot covariate data from RasterStack
#'
#' @import ggplot2 
#' 
#' @param x Object to be plotted
#' @name plot_covariate_data

plot_covariate_data <- function(x) {
  
  stopifnot(inherits(x, c('RasterStack', 'RasterBrick')))

  # cov_df <- raster::as.data.frame(x, xy = TRUE)
  # cov_df <- reshape2::melt(cov_df, id.vars = c('x','y'))
  # 
  # p <- ggplot() + geom_raster(data = cov_df, aes(x = x, y = y, fill = value)) + facet_wrap(~variable)

  p <- sp::spplot(x)
  
  print(p)
  return(p)
}

#' Plot inla.mesh object
#'
#' @param x Object to be plotted
#' @name plot_inla_mesh

plot_inla_mesh <- function(x) {
  
  stopifnot(inherits(x, 'inla.mesh'))
  
  p <- INLA::plot.inla.mesh(x)
  
  print(p)
  return(p)
  
}