#' Summary function for disaggregation fit result
#'
#' Function that summarises the result of the fit from the disaggregation model.
#'
#' Prints the negative log likelihood, model parameters and calculates metrics from in-sample performance.
#'
#' @param object Object returned from disag_model.
#' @param ... Further arguments to \emph{summary} function.
#'
#' @return A list of the model parameters, negative log likelihood and metrics from in-sample performance.
#'
#' @method summary disag_model
#'
#' @export
#' @importFrom stats cor quantile sd

summary.disag_model <- function(object, ...) {

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

  cat(paste('Likelihood function:', object$model_setup$family, '\n'))
  cat(paste('Link function:', object$model_setup$link, '\n'))

  cat('Model parameters:\n')
  print(model_params)

  cat(paste0('\nModel convergence: ', object$opt$convergence, ' (', object$opt$message, ')'))

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
#' Function that prints the result of the fit from the disaggregation model.
#'
#' Prints the negative log likelihood, model parameters and calculates metrics from in-sample performance.
#'
#' @param x Object returned from disag_model.
#' @param ... Further arguments to \emph{print} function.
#'
#' @return NULL
#'
#' @method print disag_model
#'
#' @export
#' @importFrom stats cor quantile sd

print.disag_model <- function(x, ...){

  model_params <- summary(x$sd_out, select = 'fixed')

  cat('Bayesian disaggregation model result\n')
  cat('\n')
  cat(paste('Likelihood function:', x$model_setup$family, '\n'))
  cat(paste('Link function:', x$model_setup$link, '\n'))

  cat('\nParameter values:\n')
  print(model_params[ , 1])

  return(invisible(x))
}


#' Summary function for disaggregation input data
#'
#' Function that summarizes the input data from the disaggregation model.
#'
#' Prints the number of polyons and pixels, the number of pixels in the largest and smallest polygons and summaries of the covariates.
#'
#' @param object Object returned from prepare_data.
#' @param ... Further arguments to \emph{summary} function.
#'
#' @return A list of the number of polyons, the number of covariates and summaries of the covariates.
#'
#' @method summary disag_data
#'
#' @export

summary.disag_data <- function(object, ...) {

  n_polygons <- nrow(object$polygon_shapefile)
  n_covariates <- as.integer(terra::nlyr(object$covariate_rasters))

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
#' @param x Object returned from prepare_data.
#' @param ... Further arguments to \emph{print} function.
#'
#' @return NULL
#'
#' @method print disag_data
#'
#' @export

print.disag_data <- function(x, ...){

  n_polygons <- nrow(x$polygon_shapefile)
  n_covariates <- terra::nlyr(x$covariate_rasters)

  cat(paste("They data contains", n_polygons, "polygons and", nrow(x$covariate_data), "pixels\n"))

  cat(paste("The largest polygon contains", max(table(x$covariate_data[ , x$shapefile_names$id_var])), "pixels",
            "and the smallest polygon contains", min(table(x$covariate_data[ , x$shapefile_names$id_var])), "pixels\n"))

  cat(paste("There are", n_covariates, "covariates\n"))

  return(invisible(x))
}


#' Summary function for disaggregation prediction
#'
#' Function that summarizes the prediction from the disaggregation model.
#'
#' Prints the number of polyons and pixels, the number of pixels in the largest and smallest polygons and summaries of the covariates.
#'
#' @param object Object returned from predict.disag_model
#' @param ... Further arguments to \emph{summary} function.
#'
#' @return A list of the number of polyons, the number of covariates and summaries of the covariates.
#'
#' @method summary disag_prediction
#'
#' @export

summary.disag_prediction <- function(object, ...) {

  number_realisations <- as.integer(terra::nlyr(object$uncertainty_prediction$realisations))
  max_mean <- max(terra::values(object$mean_prediction$prediction))
  min_mean <- min(terra::values(object$mean_prediction$prediction))
  max_iqr <- max((terra::values(object$uncertainty_prediction$predictions_ci[[2]]) - terra::values(object$uncertainty_prediction$predictions_ci[[1]])))
  min_iqr <- min((terra::values(object$uncertainty_prediction$predictions_ci[[2]]) - terra::values(object$uncertainty_prediction$predictions_ci[[1]])))

  cat('Predction from disaggregation model\n')
  cat('\n')
  cat('Components of the model: ')
  if(!is.null(object$mean_prediction$covariates)) cat('covariates ')
  if(!is.null(object$mean_prediction$field)) cat('field ')
  if(!is.null(object$mean_prediction$iid)) cat('iid ')
  cat('\n\n')
  cat(paste0('There are ', number_realisations, ' uncertainty realisations\n'))
  cat('\n')
  cat(paste('The mean predicted values range from', signif(min_mean, 3), 'to', signif(max_mean, 3), '\n'))
  cat(paste('The predicted IQR takes values from', signif(min_iqr, 3), 'to', signif(max_iqr, 3), '\n'))

  summary <- list(number_realisations = number_realisations,
                  range_mean_values = c(min_mean, max_mean),
                  range_iqr_values = c(min_iqr, max_iqr))

  return(invisible(summary))

}


#' Print function for disaggregation prediction
#'
#' Function that prints the prediction from the disaggregation model.
#'
#' Prints the number of polyons and pixels, the number of pixels in the largest and smallest polygons and summaries of the covariates.
#'
#' @param x Object returned from predict.disag_model.
#' @param ... Further arguments to \emph{print} function.
#'
#' @return NULL
#'
#' @method print disag_prediction
#'
#' @export

print.disag_prediction <- function(x, ...){

  cat('Predction from disaggregation model\n')
  cat('\n')
  cat('Components of the model: ')
  if(!is.null(x$mean_prediction$covariates)) cat('covariates ')
  if(!is.null(x$mean_prediction$field)) cat('field ')
  if(!is.null(x$mean_prediction$iid)) cat('iid ')
  cat('\n\n')
  cat(paste0('There are ', terra::nlyr(x$uncertainty_prediction$realisations), ' uncertainty realisations'))

  return(invisible(x))
}
