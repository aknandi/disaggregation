#' Function to predict mean from the model result
#' 
#' @param model_output fit.result object returned by fit_model function
#' 
#' @name predict_model
#'
#' @examples 
#' \dontrun{
#' predict_model(data, result)
#' }
#' 
#' @export

predict_model <- function(model_output) {
  
  data <- model_output$data

  coords <- getCoords(data)
  Amatrix <- getAmatrix(data$mesh, coords)
  
  pars <- model_output$obj$env$last.par.best
  pars <- split(pars, names(pars))
  
  # Extract field values
  field <- (Amatrix %*% pars$nodemean)[, 1]
  field_ras <- raster::rasterFromXYZ(cbind(coords, field))
  
  # Create linear predictor
  covs_by_betas <- list()
  for(i in seq_len(raster::nlayers(data$covariate_rasters))){
    covs_by_betas[[i]] <- pars$slope[i] * data$covariate_rasters[[i]]
  }
  
  cov_by_betas <- raster::stack(covs_by_betas)
  cov_contribution <- sum(cov_by_betas) + pars$intercept
  
  linear_pred <- cov_contribution + field_ras
  
  mean_prediction <- 1 / (1 + exp(-1 * linear_pred))
  
  predictions <- list(prediction = mean_prediction, 
                      field = field_ras,
                      covariates = cov_contribution)
  
  class(predictions) <- c('predictions', 'list')
  
  return(predictions)
  
}

#' Function to predict uncertainty from the model result
#' 
#' @param data disag.data object returned by prepare_data function
#' @param model_output fit.result object returned by fit_model function
#' @param N number of realisations. Default: 100
#' @param CI confidence interval. Default: 0.95
#' 
#' @name predict_uncertainty
#'
#' @examples 
#' \dontrun{
#' predict_uncertainty(data, result)
#' }
#' 
#' @export

predict_uncertainty <- function(model_output, N = 100, CI = 0.95) {
  
  data <- model_output$data
  parameters <- model_output$obj$env$last.par.best
  
  ch <- Matrix::Cholesky(model_output$sd_out$jointPrecision)
  par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = TRUE)
  
  coords <- getCoords(data)
  Amatrix <- getAmatrix(data$mesh, coords)

  predictions <- list()
  
  for(r in seq_len(N)) {
    
    p <- split(par_draws[r, ], names(parameters))
    
    # Extract field values
    field <- (Amatrix %*% p$nodemean)[, 1]
    field_ras <- raster::rasterFromXYZ(cbind(coords, field))
    
    # Create linear predictor
    covs_by_betas <- list()
    for(i in seq_len(raster::nlayers(data$covariate_rasters))){
      covs_by_betas[[i]] <- p$slope[i] * data$covariate_rasters[[i]]
    }
    cov_by_betas <- raster::stack(covs_by_betas)
    cov_contribution <- sum(cov_by_betas) + p$intercept
    
    linear_pred <- cov_contribution + field_ras
    
    predictions[[r]] <- 1 / (1 + exp(-1 * linear_pred))
  }

  predictions <- raster::stack(predictions)
  
  probs <- c((1 - CI) / 2, 1 - (1 - CI) / 2)
  predictions_ci <- raster::calc(predictions, function(x) quantile(x, probs = probs, na.rm = TRUE))
  
  uncertainty <- list(realisations = predictions,
                      predictions_ci = predictions_ci)
  
  class(uncertainty) <- c('uncertainty', 'list')
  
  return(uncertainty)
}

#' Get coordinates from raster
#'
#' @param data disag.data object 
#' @name getCoords

getCoords <- function(data) {
  
  points_raster <- data$covariate_rasters[[1]]
  points_raster[is.na(points_raster)] <- -9999
  raster_pts <- raster::rasterToPoints(points_raster, spatial = TRUE)
  coords <- raster_pts@coords
  
  return(coords)
}

#' Get Amatrix for field
#'
#' @param mesh mesh used in the model fitting
#' @param coords coordinates extracted from raster
#' 
#' @name getAmatrix

getAmatrix <- function(mesh, coords) {
  
  spde <- (INLA::inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  n_s <- nrow(spde$M0)						
  
  Amatrix <- INLA::inla.mesh.project(mesh, loc = as.matrix(coords))$A
  
  return(Amatrix)
}
