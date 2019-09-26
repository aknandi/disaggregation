#' Function to predict mean from the model result
#' 
#' @param model_output fit.result object returned by fit_model function
#' @param newdata If NULL, predictions are made using the data in model_output. 
#'   If this is a raster stack or brick, predictions will be made over this data. 
#' @param predict_iid If TRUE, any polygon iid effect from the model will be used in the prediction. Default FALSE
#'
#' @name predict_model
#'
#' @examples 
#' \dontrun{
#' predict_model(result)
#' }
#' 
#' @export

predict_model <- function(model_output, newdata = NULL, predict_iid = FALSE) {
  
  newdata <- check_newdata(newdata, model_output)

  # Pull out original data
  data <- model_output$data

  # Decide which covariates to predict over
  if(is.null(newdata)){
    covariates <- data$covariate_rasters
  } else {
    covariates <- newdata
  }
  
  # If there is no iid effect in the model, it cannot be predicted
  if(!model_output$model_setup$iid) {
    predict_iid <- FALSE
  }

  data$covariate_rasters <- covariates
  coords <- getCoords(data)
  Amatrix <- getAmatrix(data$mesh, coords)
  
  pars <- model_output$obj$env$last.par.best
  pars <- split(pars, names(pars))
  
  # Create linear predictor
  covs_by_betas <- list()
  for(i in seq_len(raster::nlayers(covariates))){
    covs_by_betas[[i]] <- pars$slope[i] * covariates[[i]]
  }
  
  cov_by_betas <- raster::stack(covs_by_betas)
  cov_contribution <- sum(cov_by_betas) + pars$intercept
  linear_pred <- cov_contribution  


  if(model_output$model_setup$field){
    # Extract field values
    field <- (Amatrix %*% pars$nodemean)[, 1]
    field_ras <- raster::rasterFromXYZ(cbind(coords, field))
    linear_pred <- linear_pred + field_ras
  } else {
    field_ras <- NULL
  }
  
  if(predict_iid) {
    tmp_shp <- model_output$data$polygon_shapefile
    tmp_shp@data <- data.frame(area_id = factor(model_output$data$polygon_data$area_id))
    shapefile_raster <- raster::rasterize(tmp_shp, 
                                          model_output$data$covariate_rasters, 
                                          field = 'area_id')
    shapefile_ids <- raster::unique(shapefile_raster)
    
    iid_ras <- shapefile_raster
    for(i in seq_along(iid_samples)) {
      iid_ras@data@values[which(shapefile_raster@data@values == shapefile_ids[i])] <- pars$iideffect[i]
    }
    linear_pred <- linear_pred + iid_ras
  } else {
    iid_ras <- NULL
  }
  
  if(model_output$model_setup$link == 0) {
    mean_prediction <- 1 / (1 + exp(-1 * linear_pred))
  } else if(model_output$model_setup$link == 1) {
    mean_prediction <- exp(linear_pred)
  } else if(model_output$model_setup$link == 2) {
    mean_prediction <- linear_pred
  }
  
  predictions <- list(prediction = mean_prediction, 
                      field = field_ras,
                      iid = iid_ras,
                      covariates = cov_contribution)
  
  class(predictions) <- c('predictions', 'list')
  
  return(predictions)
  
}

#' Function to predict uncertainty from the model result
#' 
#' @param model_output fit.result object returned by fit_model function
#' @param newdata If NULL, predictions are made using the data in model_output. 
#'   If this is a raster stack or brick, predictions will be made over this data. 
#' @param predict_iid If TRUE, any polygon iid effect from the model will be used in the prediction. Default FALSE
#' @param N number of realisations. Default: 100
#' @param CI confidence interval. Default: 0.95
#' 
#' @name predict_uncertainty
#'
#' @examples 
#' \dontrun{
#' predict_uncertainty(result)
#' }
#' 
#' @export

predict_uncertainty <- function(model_output, newdata = NULL, predict_iid = FALSE, N = 100, CI = 0.95) {
  
  
  newdata <- check_newdata(newdata, model_output)

  # Pull out original data
  data <- model_output$data

  # Decide which covariates to predict over
  if(is.null(newdata)){
    covariates <- data$covariate_rasters
  } else {
    covariates <- newdata
  }

  data$covariate_rasters <- covariates
  parameters <- model_output$obj$env$last.par.best
  
  ch <- Matrix::Cholesky(model_output$sd_out$jointPrecision)
  par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = TRUE)
  
  coords <- getCoords(data)
  Amatrix <- getAmatrix(data$mesh, coords)

  if(predict_iid) {
    tmp_shp <- model_output$data$polygon_shapefile
    tmp_shp@data <- data.frame(area_id = factor(model_output$data$polygon_data$area_id))
    shapefile_raster <- raster::rasterize(tmp_shp, 
                                          model_output$data$covariate_rasters, 
                                          field = 'area_id')
    shapefile_ids <- raster::unique(shapefile_raster)
  }
  
  predictions <- list()
  
  for(r in seq_len(N)) {
    
    p <- split(par_draws[r, ], names(parameters))
    
    # Create linear predictor
    covs_by_betas <- list()
    for(i in seq_len(raster::nlayers(covariates))){
      covs_by_betas[[i]] <- p$slope[i] * covariates[[i]]
    }
    cov_by_betas <- raster::stack(covs_by_betas)
    cov_contribution <- sum(cov_by_betas) + p$intercept
    
    linear_pred <- cov_contribution  
    
    
    if(model_output$model_setup$field){
      # Extract field values
      field <- (Amatrix %*% p$nodemean)[, 1]
      field_ras <- raster::rasterFromXYZ(cbind(coords, field))
      linear_pred <- linear_pred + field_ras
    }
    
    if(predict_iid) {
      iid_samples <- p$iideffect
      iid_ras <- shapefile_raster
      for(i in seq_along(iid_samples)) {
        iid_ras@data@values[which(shapefile_raster@data@values == shapefile_ids[i])] <- iid_samples[i]
      }
      linear_pred <- linear_pred + iid_ras
    }
    
    if(model_output$model_setup$link == 0) {
      predictions[[r]] <- 1 / (1 + exp(-1 * linear_pred))
    } else if(model_output$model_setup$link == 1) {
      predictions[[r]] <- exp(linear_pred)
    } else if(model_output$model_setup$link == 2) {
      predictions[[r]] <- linear_pred
    }
  }

  predictions <- raster::stack(predictions)
  
  probs <- c((1 - CI) / 2, 1 - (1 - CI) / 2)
  predictions_ci <- raster::calc(predictions, function(x) stats::quantile(x, probs = probs, na.rm = TRUE))
  
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



# Helper to check and sort out new raster data.
check_newdata <- function(newdata, model_output){
  if(is.null(newdata)) return(NULL)
  if(!is.null(newdata)){
    if(!(inherits(newdata, c('RasterStack', 'RasterBrick', 'RasterLayer')))){
      stop('newdata should be NULL or a RasterStack or a RasterBrick')
    } 
    if(!all(names(model_output$data$covariate_rasters) %in% names(newdata))){
      stop('All covariates used to fit the model must be in newdata')
    }
    # Take just the covariates we need and in the right order
    newdata <- newdata[[names(model_output$data$covariate_rasters)]]
  }
  return(newdata)
}


