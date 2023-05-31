#' Plot input data for disaggregation
#'
#' Plotting function for class \emph{disag_data} (the input data for disaggregation).
#' 
#' Produces three plots: polygon response data, covariate rasters and INLA mesh.
#'
#' @param x Object of class \emph{disag_data} to be plotted.
#' @param which If a subset of plots is required, specify a subset of the numbers 1:3
#' @param ... Further arguments to \emph{plot} function.
#' 
#' @return A list of three plots: the polygon plot (ggplot), covariate plot (spplot) and INLA mesh plot (ggplot)
#' 
#' @import ggplot2
#' @method plot disag_data
#' 
#' @export

plot.disag_data <- function(x, which = c(1,2,3), ...) {
  
  plots <- list()
  titles <- c()
  
  if(1 %in% which) {
    plots$polygon <- plot_polygon_data(x$polygon_shapefile, x$shapefile_names)
    titles <- c(titles, 'Polygon response data')
  }
  
  if(2 %in% which) {
    stopifnot(inherits(x$covariate_rasters, c('RasterStack', 'RasterBrick')))
    plots$covariates <- sp::spplot(x$covariate_rasters)
    titles <- c(titles, 'Covariate rasters')
  }
  
  if(3 %in% which & !is.null(x$mesh)) {
    stopifnot(inherits(x$mesh, 'inla.mesh'))
    plots$mesh <- plot_mesh(x$mesh)
    titles <- c(titles, 'INLA mesh for spatial field')
  }
  
  print(cowplot::plot_grid(plotlist = plots, labels = titles, label_size = 10))
  
  return(invisible(plots))
}

#' Plot results of fitted model
#'
#' Plotting function for class \emph{disag_model} (the result of the disaggregation fitting).
#' 
#' Produces two plots: results of the fixed effects and in-sample observed vs predicted plot.
#' 
#' @param x Object of class \emph{disag_model} to be plotted.
#' @param ... Further arguments to \emph{plot} function.
#' 
#' @return A list of two ggplot plots: results of the fixed effects and an in-sample observed vs predicted plot 
#' 
#' @import ggplot2
#' @method plot disag_model
#' 
#' @export


plot.disag_model <- function(x, ...){
  
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

#' Plot mean and uncertainty predictions from the disaggregation model results
#'
#' Plotting function for class \emph{disag_prediction} (the mean and uncertainty predictions of the disaggregation fitting).
#' 
#' Produces raster plots of the mean prediction, and the lower and upper confidence intervals.
#'
#' @param x Object of class \emph{disag_prediction} to be plotted.
#' @param ... Further arguments to \emph{plot} function.
#' 
#' @return A list of plots of rasters from the prediction: mean prediction, lower CI and upper CI.
#' 
#' @method plot disag_prediction
#' 
#' @export


plot.disag_prediction <- function(x, ...) {

  rasters_to_plot <- raster::stack(x$mean_prediction$prediction, x$uncertainty_prediction$predictions_ci)
  names(rasters_to_plot) <- c('mean prediction', 'lower CI', 'upper CI')
  
  plots <- sp::spplot(rasters_to_plot)
  
  print(plots)
  
  return(invisible(plots))
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
    scale_fill_viridis_c(trans = 'identity')
  
  return(invisible(p))
}

# A ggplot2 method for plotting INLA mesh objects.
#
# @param object An inla.mesh object
# @param col Colour for data points
# @param lwd Line width
# @param linecol The colour for the mesh edges
# @param size size Size of data points
# 
# @name plot_mesh

plot_mesh <- function(x, main = '', col = 'blue', lwd = 0.5, linecol = 'darkgrey', size = 1.2) {
  
  mesh <- x
  # extract point data
  d <- data.frame(x = mesh$loc[, 1], y = mesh$loc[, 2], type = 'evertices')
  levels(d$type) <- c('evertices', 'adata')
  d[mesh$idx$loc, 'type'] <- 'adata'
  # extract lines data. 
  # mesh$graph$tv column 1, 2, 3 are points in triangles.
  # Therefore need 1 to 2, 2 to 3 and 3 to 1.
  idx = rbind(mesh$graph$tv[, 1:2, drop = FALSE], 
              mesh$graph$tv[, 2:3, drop = FALSE], 
              mesh$graph$tv[, c(3, 1), drop = FALSE])
  segments <- data.frame(mesh$loc[idx[, 1], 1:2], mesh$loc[idx[, 2], 1:2], type = 'bsegments')
  
  innerouter <- data.frame(mesh$loc[mesh$segm$bnd$idx[, 1], 1:2],
                           mesh$loc[mesh$segm$bnd$idx[, 2], 1:2],
                           type = 'cbinding', stringsAsFactors = FALSE)
  if(nrow(mesh$segm$int$idx) > 0){
    innerouter <- rbind(innerouter,
                        data.frame(mesh$loc[mesh$segm$int$idx[, 1], 1:2],
                                   mesh$loc[mesh$segm$int$idx[, 2], 1:2],
                                   type = 'dinternal'))
  } else {
    #innerouter <- rbind(innerouter,
    #                    NA)
    #innerouter[nrow(innerouter), 5] <- 'dinternal'
    innerouter$type = factor(innerouter$type, levels = c('dinternal', 'cbinding'))
  }
  
  
  names(segments) <- c('x1', 'y1', 'x2', 'y2', 'type')
  names(innerouter) <- c('x1', 'y1', 'x2', 'y2', 'type')
  
  segments <- rbind(segments, innerouter)
  
  
  p <- ggplot2::ggplot(data = d, 
                       ggplot2::aes_string('x', 'y', 
                                           colour = 'type', 
                                           size = 'type')) +
    ggplot2::geom_segment(data = segments, 
                          ggplot2::aes_string(x = 'x1', y = 'y1', xend = 'x2', yend = 'y2')) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = 'none')
  #stroke
  p <- p +
    ggplot2::scale_colour_manual(values = c(col, linecol, 'black', 'black', 'black'), drop = FALSE) +
    ggplot2::scale_size_manual(values = c(size, lwd, 1.3, 1.3, 0), drop = FALSE) +
    ggtitle(main)
  
  return(invisible(p))
}
