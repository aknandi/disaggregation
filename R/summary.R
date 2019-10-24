#' Summary function for disaggregation fit result
#' 
#' Function that summarises the fit result from the disaggregation model.
#' 
#' Prints the negative log likelihood, model parameters and calculates metrics from in-sample performance.
#'
#' @param object Object returned from fit_model.
#' @param ... Further arguments to \emph{summary} function.
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
  if(object$model_setup$family == 'gaussian') {
    observed_data = report$polygon_response_data/report$reportnormalisation
    predicted_data = report$reportprediction_rate
  } else if(object$model_setup$family == 'binomial') {
    observed_data = object$data$polygon_data$response/object$data$polygon_data$N
    predicted_data = report$reportprediction_rate
  } else if(object$model_setup$family == 'poisson') {
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


#' Print function for disaggregation fit result.
#' 
#' Function that prints the fit result from the disaggregation model.
#' 
#' Prints the negative log likelihood, model parameters and calculates metrics from in-sample performance.
#'
#' @param x Object returned from fit_model.
#' @param ... Further arguments to \emph{print} function.
#' 
#' @method print fit.result
#' 
#' @export
#' @importFrom stats cor quantile sd

print.fit.result <- function(x, ...){
  summary(x)
  return(NULL)
  
}


#' Summary function for disaggregation input data
#' 
#' Function that summarizes the input data from the disaggregation model.
#' 
#' Prints the number of polyons and pixels, the number of pixels in the largest and smallest polygons and summaries of the covariates.
#'
#' @param object Object returned from fit_model.
#' @param ... Further arguments to \emph{summary} function.
#' 
#' @method summary disag.data
#' 
#' @export

summary.disag.data <- function(object, ...) {

  n_polygons <- nrow(object$polygon_shapefile)
  n_covariates <- raster::nlayers(object$covariate_rasters)
  covariate_names <- names(object$covariate_rasters)
  
  cat(paste("They data contains", n_polygons, "polygons and", nrow(object$covariate_data), "pixels\n"))
  
  cat(paste("The largest polygon contains", max(table(object$covariate_data[ , object$shapefile_names$id_var])), "pixels", 
            "and the smallest polygon contains", min(table(object$covariate_data[ , object$shapefile_names$id_var])), "pixels\n"))
  
  cat(paste("There are", n_covariates, "covariates\n"))
  
  covariate_summary <- summary(object$covariate_data[ , names(object$covariate_rasters)])
  
  cat("\nCovariate summary:\n")
  print(covariate_summary)
  
  summary <- list(number_polygons = n_polygons,
                  number_covariates = n_covariates,
                  covariate_summary = covariate_summary)
  
  return(invisible(summary))
  
}


#' Print function for disaggregation input data
#' 
#' Function that prints the input data from the disaggregation model.
#' 
#' Prints the number of polyons and pixels, the number of pixels in the largest and smallest polygons and summaries of the covariates.
#'
#' @param x Object returned from fit_model.
#' @param ... Further arguments to \emph{print} function.
#' 
#' @method print disag.data
#' 
#' @export
print.disag.data <- function(x, ...){
  summary(x)
  return(NULL)
}
