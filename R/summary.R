#' Summary function for fit result
#'
#' @param object Object returned from fit_model
#' @param ... Further arguments passed to or from other methods.
#' 
#' @method summary fit.result
#' 
#' @export

summary.fit.result <- function(object, ...) {
  
  model_params <- summary(object$sd_out, select = 'fixed')
  
  report <- object$obj$report()
  nll <- report$nll
  
  in_sample <- data.frame(obs = report$polygon_response_data, pred = report$reportprediction)
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
  
  return(summary)
  
}