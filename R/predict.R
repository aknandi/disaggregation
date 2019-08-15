#' Function to predict from the model result
#' 
#' @param data disag.data object returned by prepare_data function
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

predict_model <- function(data, model_output) {
  
  pars <- model_output$obj$env$last.par.best
  
  # Extract the Amatrix and coords of the field
  points_raster <- data$covariate_rasters[[1]]
  points_raster[is.na(points_raster)] <- -9999
  raster_pts <- raster::rasterToPoints(points_raster, spatial = TRUE)
  coords <- raster_pts@coords
  
  # Get random field predicted
  spde <- (INLA::inla.spde2.matern(data$mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  n_s <- nrow(spde$M0)						
  
  Amatrix <- INLA::inla.mesh.project(data$mesh, loc = as.matrix(coords))$A
  
  # Split up parameters
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