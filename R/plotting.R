#' Plot input data for disaggregation
#'
#' Plotting function for class \emph{disag.data} (the input data for disaggragation).
#' 
#' Produces three plots: polygon response data, covariate rasters and INLA mesh.
#'
#' @param x Object of class \emph{disag.data} to be plotted.
#' @param ... Further arguments to \emph{plot} function.
#' 
#' @return A list of two plots: the polygon plot (ggplot) and covariate plot (spplot)
#' 
#' @import ggplot2
#' @method plot disag.data
#' 
#' @export

plot.disag.data <- function(x, ...) {
  
  plots <- list()
  
  plots$polygon <- plot_polygon_data(x$polygon_shapefile, x$shapefile_names)
  
  stopifnot(inherits(x$covariate_rasters, c('RasterStack', 'RasterBrick')))
  plots$covariates <- sp::spplot(x$covariate_rasters, main = 'Covariate rasters')
  print(plots$covariates)
  
  if(!is.null(x$mesh)) {
    stopifnot(inherits(x$mesh, 'inla.mesh'))
    INLA::plot.inla.mesh(x$mesh, main = 'INLA mesh for spatial field')
  }
  
  return(invisible(plots))
}

#' Plot results of fitted model
#'
#' Plotting function for class \emph{fit.result} (the result of the disaggragation fitting).
#' 
#' Produces two plots: results of the fixed effects and in-sample observed vs predicted plot.
#' 
#' @param x Object of class \emph{fit.result} to be plotted.
#' @param ... Further arguments to \emph{plot} function.
#' 
#' @return A list of two ggplot plots: results of the fixed effects and an in-sample observed vs predicted plot 
#' 
#' @import ggplot2
#' @method plot fit.result
#' 
#' @export


plot.fit.result <- function(x, ...){
  
  parameter <- sd <- obs <- pred <- NULL
  posteriors <- as.data.frame(summary(x$sd_out, select = 'fixed'))
  posteriors <- dplyr::mutate(posteriors, name = rownames(posteriors))
  names(posteriors) <- c('mean', 'sd', 'parameter')

  fixedeffects <- ggplot() + 
    geom_errorbar(posteriors, mapping = aes(x = parameter, ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "blue") + 
    geom_point(posteriors, mapping = aes(x = parameter, y = mean)) + 
    ggtitle("Fixed effects")
  
  report <- x$obj$report()
  
  # Form of the observed and predicted results depends on the likelihood function used
  if(x$model_setup$family == 'gaussian') {
    observed_data = report$polygon_response_data/report$reportnormalisation
    predicted_data = report$reportprediction_rate
    title <- 'In sample performance: incidence rate'
  } else if(x$model_setup$family == 'binomial') {
    observed_data = x$data$polygon_data$response/x$data$polygon_data$N
    predicted_data = report$reportprediction_rate
    title <- 'In sample performance: prevalence rate'
  } else if(x$model_setup$family == 'poisson') {
    observed_data = report$polygon_response_data/report$reportnormalisation
    predicted_data = report$reportprediction_rate
    title <- 'In sample performance: incidence rate'
  }
  
  data <- data.frame(obs = observed_data, pred = predicted_data)
  
  obspred <- ggplot(data, aes(x = obs, y = pred)) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, color = 'blue') + 
    ggtitle(title)
  
  plots <- list(fixedeffects, obspred)
  print(cowplot::plot_grid(plotlist = plots))
  
  return(invisible(plots))
}

#' Plot predictions from the disaggregation model results
#'
#' Plotting function for class \emph{predictions} (the mean predictions of the disaggragation fitting).
#' 
#' Produces plots of the mean prediction, and the covariate, field  and iid contribution to the linear predictor.
#'
#' @param x Object of class \emph{predictions} to be plotted.
#' @param ... Further arguments to \emph{plot} function.
#' 
#' @return A list of plots of rasters from the prediction: mean prediction, covariate contribution, field contribution (if used) 
#' and the iid contribution (if used)
#' 
#' @method plot predictions
#' 
#' @export


plot.predictions <- function(x, ...) {

  mean_plot <- sp::spplot(x$prediction, main = list(label = 'mean prediction'))
  print(mean_plot)
  
  covariate_plot <- sp::spplot(x$covariates, main = list(label = 'covariate contribution'))
  print(covariate_plot)
  
  plots <- list(mean = mean_plot, covariates = covariate_plot)
  
  if(!is.null(x$field)) {
    field_plot <- sp::spplot(x$field, main = list(label = 'spatial field'))
    plots <- c(plots, field = list(field_plot))
    print(field_plot)
  }
  if(!is.null(x$iid)) {
    iid_plot <- sp::spplot(x$iid, main = list(label = 'iideffect'))
    plots <- c(plots, iid = list(iid_plot))
    print(iid_plot)
  }
  
  return(invisible(plots))
}

#' Plot uncertainty predictions from the disaggregation model results
#'
#' Plotting function for class \emph{uncertainty} (the uncertainty predictions of the disaggragation fitting).
#' 
#' Produces a plot of the lower and upper credible interval rasters.
#' 
#' @param x Object of class \emph{uncertainty} to be plotted.
#' @param ... Further arguments to \emph{plot} function.
#' 
#' @return A plot of the lower and upper credible interval rasters
#' 
#' @method plot uncertainty
#' 
#' @export


plot.uncertainty <- function(x, ...) {
  
  unc_plot <- sp::spplot(x$predictions_ci, main = 'uncertainty predictions')
  
  print(unc_plot)
  
  return(invisible(unc_plot))
}

# Plot polygon data from SpatialPolygonDataFrame
#
# @param x Object to be plotted
# @param names list of 2 names: polygon id variable and response variable names
# 
# @return A ggplot of the polygon data
# 
# @name plot_polygon_data

plot_polygon_data <- function(x, names) {

  # Rename the response variable for plotting
  shp <- x
  shp@data <- dplyr::rename(shp@data, 'response' = names$response_var)
  shp@data <- dplyr::rename(shp@data, 'area_id' = names$id_var)
  
  area_id <- long <- lat <- group <- response <- NULL
  stopifnot(inherits(shp, 'SpatialPolygonsDataFrame'))
  
  df_fortify <- fortify(shp, region = 'area_id')
  
  df <- shp@data
  df <- dplyr::mutate(df, area_id = as.character(area_id)) 
  df <- dplyr::left_join(df_fortify, df, by = c('id' = 'area_id'))
  
  p <- ggplot(df, aes(long, lat, group = group, fill = response)) + 
    geom_polygon() +
    coord_equal() +
    scale_fill_viridis_c(trans = 'identity') +
    ggtitle('Polygon response data')
  
  print(p)
  return(p)
}


