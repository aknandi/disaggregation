#' These are plot methods for input data and predictions for disaggragation
#'
#' @param x Object to be plotted
#' @param ... Further arguments passed to or from other methods.
#' 
#' @import ggplot2
#' @method plot disag.data
#' 
#' @export

plot.disag.data <- function(x, ...) {
  
  # Plot polygon data, covariate rasters and mesh
  plots <- list()
  
  plots$polygon <- plot_polygon_data(x$polygon_shapefile)
  plots$covariates <- plot_covariate_data(x$covariate_rasters)
  plots$mesh <- plot_inla_mesh(x$mesh)
  
  return(invisible(plots))
}

#' Plot results of fitted model
#'
#' @param x Object to be plotted
#' @param ... Further arguments passed to or from other methods.
#' 
#' @import ggplot2
#' @method plot fit.result
#' 
#' @export


plot.fit.result <- function(x, ...){
  
  posteriors <- as.data.frame(summary(x$sd_out, select = 'fixed'))
  posteriors <- dplyr::mutate(posteriors, name = rownames(posteriors))
  names(posteriors) <- c('mean', 'sd', 'parameter')

  fixedeffects <- ggplot() + 
    geom_errorbar(posteriors, mapping = aes(x = parameter, ymin = mean - sd, ymax = mean + sd), width=0.2, color="blue") + 
    geom_point(posteriors, mapping = aes(x = parameter, y = mean)) + 
    ggtitle("Fixed effects")
  
  report <- x$obj$report()
  data <- data.frame(obs = report$polygon_coverage_data, pred = report$reportprediction)
  
  obspred <- ggplot(data, aes(x = obs, y = pred)) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, color = 'blue') + 
    ggtitle("In sample performance")
  
  plots <- list(fixedeffects, obspred)
  print(cowplot::plot_grid(plotlist = plots))
  
  return(invisible(plots))
}

#' Plot predictions
#'
#' @param x Object to be plotted
#' @param ... Further arguments passed to or from other methods.
#' 
#' @method plot predictions
#' 
#' @export


plot.predictions <- function(x, ...) {
  
  mean_plot <- sp::spplot(x$prediction)
  field_plot <- sp::spplot(x$field)
  covariate_plot <- sp::spplot(x$covariates)
  
  print(mean_plot)
  print(field_plot)
  print(covariate_plot)
  
  plots <- list(mean_plot, field_plot, covariate_plot)
  
  return(invisible(plots))
}

#' Plot uncertainty
#' 
#' @param x Object to be plotted
#' @param ... Further arguments passed to or from other methods.
#' 
#' @method plot uncertainty
#' 
#' @export


plot.uncertainty <- function(x, ...) {
  
  unc_plot <- sp::spplot(x$predictions_ci)
  
  print(unc_plot)
  
  return(invisible(unc_plot))
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