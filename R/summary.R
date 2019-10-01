#' Summary function for fit result
#'
#' @param object Object returned from fit_model
#' @param ... Further arguments passed to or from other methods.
#' 
#' @method summary fit.result
#' 
#' @export
#' @importFrom stats cor quantile sd

summary.fit.result <- function(object, ...) {
  
  pred <- obs <- NULL
  
  model_params <- summary(object$sd_out, select = 'fixed')
  
  report <- object$obj$report()
  nll <- report$nll
  
  # Form of the observed and predicted results depends on the likelihood function used
  if(object$model_setup$family == 0) {
    observed_data = report$polygon_response_data/report$reportnormalisation
    predicted_data = report$reportprediction_rate
  } else if(object$model_setup$family == 1) {
    observed_data = object$data$polygon_data$response/object$data$polygon_data$N
    predicted_data = report$reportprediction_rate
  } else if(object$model_setup$family == 2) {
    observed_data = report$polygon_response_data
    predicted_data = report$reportprediction_cases
  }
  
  in_sample <- data.frame(obs = observed_data, pred = predicted_data)
  in_sample_reduced <- in_sample[!is.na(in_sample$pred), ]
  metrics <- dplyr::summarise(in_sample_reduced, 
                              RMSE = sqrt(mean((pred - obs) ^ 2)),
                              MAE = mean(abs(pred - obs)),
                              pearson = cor(pred, obs, method = 'pearson'),
                              spearman = cor(pred, obs, method = 'spearman'),
                              log_pearson = cor(log1p(pred), log1p(obs), method = 'pearson'))
  
  cat('Model parameters:\n')
  print(model_params)
  
  cat(paste('\nNegative log likelihood: ', nll, '\n'))
  
  cat('\nIn sample performance:\n')
  print(metrics)
  
  summary <- list(model_params = model_params,
                  nll = nll,
                  metrics = metrics)
  
  return(invisible(summary))
  
}




#' Print function for fit result
#'
#' @param x Object returned from fit_model
#' @param ... Further arguments passed to or from other methods.
#' 
#' @method print fit.result
#' 
#' @export
#' @importFrom stats cor quantile sd

print.fit.result <- function(x, ...) {
  
  pred <- obs <- NULL
  
  model_params <- summary(x$sd_out, select = 'fixed')
  
  report <- x$obj$report()
  nll <- report$nll
  
  # Form of the observed and predicted results depends on the likelihood function used
  if(x$model_setup$family == 0) {
    observed_data = report$polygon_response_data/report$reportnormalisation
    predicted_data = report$reportprediction_rate
  } else if(x$model_setup$family == 1) {
    observed_data = x$data$polygon_data$response/x$data$polygon_data$N
    predicted_data = report$reportprediction_rate
  } else if(x$model_setup$family == 2) {
    observed_data = report$polygon_response_data
    predicted_data = report$reportprediction_cases
  }
  
  in_sample <- data.frame(obs = observed_data, pred = predicted_data)
  in_sample_reduced <- in_sample[!is.na(in_sample$pred), ]
  metrics <- dplyr::summarise(in_sample_reduced, 
                              RMSE = sqrt(mean((pred - obs) ^ 2)),
                              MAE = mean(abs(pred - obs)),
                              pearson = cor(pred, obs, method = 'pearson'),
                              spearman = cor(pred, obs, method = 'spearman'),
                              log_pearson = cor(log1p(pred), log1p(obs), method = 'pearson'))
  
  cat('Model parameters:\n')
  print(model_params)
  
  cat(paste('\nNegative log likelihood: ', nll, '\n'))
  
  cat('\nIn sample performance:\n')
  print(metrics)
  
  return(NULL)
  
}




#' Summary function for disag.data
#'
#' @param object Object returned from fit_model
#' @param ... Further arguments passed to or from other methods.
#' 
#' @method summary disag.data
#' 
#' @export

summary.disag.data <- function(object, ...) {

  n_polygons <- nrow(object$polygon_shapefile)
  n_covariates <- raster::nlayers(object$covariate_rasters)
  covariate_names <- names(object$covariate_rasters)
  
  cat(paste("They data contains", n_polygons, "polygons and", nrow(object$covariate_data), "pixels\n"))
  
  cat(paste("The largest polygon contains", max(table(object$covariate_data$area_id)), "pixels", 
            "and the smallest polygon contains", min(table(object$covariate_data$area_id)), "pixels\n"))
  
  cat(paste("There are", n_covariates, "covariates\n"))
  
  covariate_summary <- summary(object$covariate_data[ , names(object$covariate_rasters)])
  
  return(list(number_polygons = n_polygons,
              number_covariates = n_covariates,
              covariate_summary = covariate_summary))
  
}


#' Print method for disag.data
#'
#' @param x Object returned from fit_model
#' @param ... Further arguments passed to or from other methods.
#' 
#' @method print disag.data
#' 
#' @export

print.disag.data <- function(x, ...){
  
  n_polygons <- nrow(x$polygon_shapefile)
  n_covariates <- raster::nlayers(x$covariate_rasters)
  covariate_names <- names(x$covariate_rasters)
  
  cat(paste("They data contains", n_polygons, "polygons and", nrow(x$covariate_data), "pixels\n"))
  
  cat(paste("The largest polygon contains", max(table(x$covariate_data$area_id)), "pixels", 
            "and the smallest polygon contains", min(table(x$covariate_data$area_id)), "pixels\n"))
  
  cat(paste("There are", n_covariates, "covariates\n"))
  
  covariate_summary <- print(summary(x$covariate_data[ , names(x$covariate_rasters)]))
  
  return(NULL)
}
