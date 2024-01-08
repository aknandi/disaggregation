#' Predict mean and uncertainty from the disaggregation model result
#'
#' \emph{predict.disag_model} function takes a \emph{disag_model} object created by \emph{disaggregation::disag_model} and
#' predicts mean and uncertainty maps.
#'
#' To predict over a different spatial extent to that used in the model,
#' a SpatRaster covering the region to make predictions over is passed to the argument \emph{newdata}.
#' If this is not given predictions are made over the data used in the fit.
#'
#' The \emph{predict_iid} logical flag should be set to TRUE if the results of the iid effect from the model are to be used in the prediction.
#'
#' For the uncertainty calculations, the number of the realisations and the size of the confidence interval to be calculated
#' are given by the arguments \emph{N} and \emph{CI} respectively.
#'
#' @param object disag_model object returned by disag_model function.
#' @param newdata If NULL, predictions are made using the data in model_output.
#'   If this is a raster stack or brick, predictions will be made over this data.
#' @param predict_iid logical. If TRUE, any polygon iid effect from the model will be used in the prediction. Default FALSE.
#' @param N Number of realisations. Default: 100.
#' @param CI Confidence interval to be calculated from the realisations. Default: 0.95.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \emph{disag_prediction} which consists of a list of two objects:
#'  \item{mean_prediction }{List of:
#'   \itemize{
#'    \item \emph{prediction} Raster of mean predictions based.
#'    \item \emph{field} Raster of the field component of the linear predictor.
#'    \item \emph{iid} Raster of the iid component of the linear predictor.
#'    \item \emph{covariates} Raster of the covariate component of the linear predictor.
#'   }}
#'  \item{uncertainty_prediction: }{List of:
#'   \itemize{
#'    \item \emph{realisations} SpatRaster of realisations of predictions. Number of realisations defined by argument \emph{N}.
#'    \item \emph{predictions_ci} SpatRaster of the upper and lower credible intervals. Defined by argument \emph{CI}.
#'   }}
#'
#'
#' @method predict disag_model
#'
#' @examples
#' \dontrun{
#' predict(fit_result)
#' }
#'
#' @export


predict.disag_model <- function(object, newdata = NULL, predict_iid = FALSE, N = 100, CI = 0.95, ...) {

  mean_prediction <- predict_model(object, newdata = newdata, predict_iid)

  uncertainty_prediction <- predict_uncertainty(object, newdata = newdata, predict_iid, N, CI)

  prediction <- list(mean_prediction = mean_prediction,
                     uncertainty_prediction = uncertainty_prediction)

  class(prediction) <- c('disag_prediction', 'list')

  return(prediction)
}

#' Function to predict mean from the model result
#'
#' \emph{predict_model} function takes a \emph{disag_model} object created by
#' \emph{disaggregation::disag_model} and predicts mean maps.
#'
#' Function returns rasters of the mean predictions as well as the  covariate and field contributions
#' to the linear predictor.
#'
#' To predict over a different spatial extent to that used in the model,
#' a SpatRaster covering the region to make predictions over is passed to the argument \emph{newdata}.
#' If this is not given predictions are made over the data used in the fit.
#'
#' The \emph{predict_iid} logical flag should be set to TRUE if the results of the iid effect from the model are to be used in the prediction.
#'
#' @param model_output disag_model object returned by disag_model function
#' @param newdata If NULL, predictions are made using the data in model_output.
#'   If this is a raster stack or brick, predictions will be made over this data. Default NULL.
#' @param predict_iid If TRUE, any polygon iid effect from the model will be used in the prediction. Default FALSE.
#'
#' @return The mean prediction, which is a list of:
#'   \itemize{
#'    \item \emph{prediction} Raster of mean predictions based.
#'    \item \emph{field} Raster of the field component of the linear predictor.
#'    \item \emph{iid} Raster of the iid component of the linear predictor.
#'    \item \emph{covariates} Raster of the covariate component of the linear predictor.
#'   }
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

  objects_for_prediction <- setup_objects(model_output, newdata = newdata, predict_iid)

  pars <- model_output$obj$env$last.par.best
  pars <- split(pars, names(pars))

  prediction <- predict_single_raster(pars,
                                      objects_for_prediction,
                                      link_function = model_output$model_setup$link)

  return(prediction)

}

#' Function to predict uncertainty from the model result
#'
#' \emph{predict_uncertainty} function takes a \emph{disag_model} object created by
#' \emph{disaggregation::disag_model} and predicts upper and lower credible interval maps.
#'
#' Function returns a SpatRaster of the realisations as well as the upper and lower credible interval rasters.
#'
#' To predict over a different spatial extent to that used in the model,
#' a SpatRaster covering the region to make predictions over is passed to the argument \emph{newdata}.
#' If this is not given predictions are made over the data used in the fit.
#'
#' The \emph{predict_iid} logical flag should be set to TRUE if the results of the iid effect from the model are to be used in the prediction.
#'
#' The number of the realisations and the size of the confidence interval to be calculated.
#' are given by the arguments \emph{N} and \emph{CI} respectively.
#'
#' @param model_output disag_model object returned by disag_model function.
#' @param newdata If NULL, predictions are made using the data in model_output.
#'   If this is a raster stack or brick, predictions will be made over this data. Default NULL.
#' @param predict_iid If TRUE, any polygon iid effect from the model will be used in the prediction. Default FALSE.
#' @param N number of realisations. Default: 100.
#' @param CI confidence interval. Default: 0.95.
#'
#' @return The uncertainty prediction, which is a list of:
#'   \itemize{
#'    \item \emph{realisations} SpatRaster of realisations of predictions. Number of realisations defined by argument \emph{N}.
#'    \item \emph{predictions_ci} SpatRaster of the upper and lower credible intervals. Defined by argument \emph{CI}.
#'   }
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

  objects_for_prediction <- setup_objects(model_output, newdata = newdata, predict_iid)

  parameters <- model_output$obj$env$last.par.best

  # If we have either of the random effects, we have the jointPrecision matrix.
  #   but if we have neither, we don't get that matrix and should use the
  #   covariance matrix instead

  #CH <- Matrix::Cholesky(as(S, 'dsCMatrix'))
  #x <- rmvn.sparse(10, mu, CH, prec=FALSE) ## 10 random draws of x
  #d <- dmvn.sparse(x, mu, CH, prec=FALSE) ## densities of the 10 draws


  if(model_output$model_setup$iid | model_output$model_setup$field){
    ch <- Matrix::Cholesky(model_output$sd_out$jointPrecision)
    par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = TRUE)
  } else {
    covariance_matrix <- Matrix::Matrix(model_output$sd_out$cov.fixed, sparse = TRUE)
    ch <- Matrix::Cholesky(covariance_matrix)
    par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = FALSE)
  }

  predictions <- list()

  for(r in seq_len(N)) {

    p <- split(par_draws[r, ], names(parameters))

    prediction_result <- predict_single_raster(p,
                                               objects_for_prediction,
                                               link_function = model_output$model_setup$link)

    predictions[[r]] <- prediction_result$prediction
  }

  predictions <- terra::rast(predictions)

  probs <- c((1 - CI) / 2, 1 - (1 - CI) / 2)
  predictions_ci <- terra::app(predictions, function(x) stats::quantile(x, probs = probs, na.rm = TRUE))


  names(predictions_ci) <- c('lower CI', 'upper CI')

  uncertainty <- list(realisations = predictions,
                      predictions_ci = predictions_ci)

  return(uncertainty)
}

# Get coordinates from raster
#
# @param data disag_data object
#
# @return A matrix of the coordinates of the raster
#
# @name getCoords

getCoords <- function(data) {

  points_raster <- data$covariate_rasters[[1]]
  points_raster[is.na(terra::values(points_raster, mat = FALSE))] <- -9999
  raster_pts <- terra::as.points(points_raster)
  coords <- terra::crds(raster_pts)

    return(coords)
}

# Helper to check and sort out new raster data.
check_newdata <- function(newdata, model_output){
  if(is.null(newdata)) return(NULL)
  if(!is.null(newdata)){
    if(!(inherits(newdata, c('SpatRaster')))){
      stop('newdata should be NULL or a SpatRaster')
    }
    if(!all(names(model_output$data$covariate_rasters) %in% names(newdata))){
      stop('All covariates used to fit the model must be in newdata')
    }
    # Take just the covariates we need and in the right order
    newdata <- newdata[[names(model_output$data$covariate_rasters)]]
  }
  return(newdata)
}

# Function to setup covariates, field and iid objects for prediction
setup_objects <- function(model_output, newdata = NULL, predict_iid = FALSE) {

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

  # If there is no iid effect in the model, it cannot be predicted
  if(!model_output$model_setup$iid) {
    predict_iid <- FALSE
  }

  if(model_output$model_setup$field) {
    if(is.null(newdata)) {
      coords <- data$coordsForPrediction
    } else {
      coords <- getCoords(data)
    }
    Amatrix <- fmesher::fm_evaluator(data$mesh, loc = as.matrix(coords))$proj$A
    field_objects <- list(coords = coords, Amatrix = Amatrix)
  } else {
    field_objects <- NULL
  }

  if(predict_iid) {
    tmp_shp <- model_output$data$polygon_shapefile
    #needed to avoid errors in testing
    if (!("area_id" %in% names(model_output$data$polygon_shapefile))){
      tmp_shp <- dplyr::bind_cols(tmp_shp,
                                  area_id =
                                    factor(model_output$data$polygon_data$area_id))
    }
    shapefile_raster <- terra::rasterize(tmp_shp,
                                          model_output$data$covariate_rasters,
                                          field = 'area_id')
    shapefile_ids <- terra::unique(shapefile_raster)
    iid_objects <- list(shapefile_raster = shapefile_raster, shapefile_ids = shapefile_ids)
  } else {
    iid_objects <- NULL
  }
  return(list(covariates = covariates,
              field_objects = field_objects,
              iid_objects = iid_objects))
}

# Function to take model parameters and predict a single raster
predict_single_raster <- function(model_parameters, objects, link_function) {

  # Create linear predictor
  covs_by_betas <- list()
  for(i in seq_len(terra::nlyr(objects$covariates))){
    covs_by_betas[[i]] <- model_parameters$slope[i] * objects$covariates[[i]]
  }

  cov_by_betas <- terra::rast(covs_by_betas)
  if(terra::nlyr(cov_by_betas) > 1){
    sum_cov_by_betas <- sum(cov_by_betas)
  } else {
    # With only 1 covariate, there's nothing to sum. Do this to avoid warnings.
    sum_cov_by_betas <- cov_by_betas
  }
  cov_contribution <- sum_cov_by_betas + model_parameters$intercept

  linear_pred <- cov_contribution

  if(!is.null(objects$field_objects)){
    # Extract field values
    field <- (objects$field_objects$Amatrix %*% model_parameters$nodemean)[, 1]
    field_ras <- terra::rast(cbind(objects$field_objects$coords, field),
                             type = 'xyz',
                             crs = terra::crs(linear_pred))
    linear_pred <- linear_pred + field_ras
  } else {
    field_ras <- NULL
  }

    if(!is.null(objects$iid_objects)) {
    iid_ras <- objects$iid_objects$shapefile_raster
    iideffect_sd <- 1/sqrt(exp(model_parameters$iideffect_log_tau))
    # todo
    for(i in seq_along(model_parameters$iideffect)) {
      targetvals <- terra::values(objects$iid_objects$shapefile_raster,
                                  dataframe = FALSE, mat = FALSE)
      whichvals <- which(targetvals == objects$iid_objects$shapefile_ids[1, i])
      terra::values(iid_ras)[whichvals] <-
        model_parameters$iideffect[i]
      na_pixels <- which(is.na(terra::values(iid_ras, dataframe = FALSE, mat = FALSE)))
      na_iid_values <- stats::rnorm(length(na_pixels), 0, iideffect_sd)
      terra::values(iid_ras)[na_pixels] <- na_iid_values
    }
    if(terra::ext(iid_ras) != terra::ext(linear_pred)) {
      # Extent of prediction space is different to the original model. Keep any overlapping iid values but predict to the new extent
      raster_new_extent <- linear_pred
      terra::values(raster_new_extent) <- NA
      #iid_ras <- terra::merge(iid_ras, raster_new_extent, ext = terra::ext(raster_new_extent))
      # NOt sure why we no longer need the ext argument
      # SS - added a crop which I think does the same thing
      iid_ras <- terra::merge(iid_ras, raster_new_extent)
      iid_ras <- terra::crop(iid_ras, raster_new_extent)
      missing_pixels <- which(is.na(terra::values(iid_ras, dataframe = FALSE, mat = FALSE)))
      missing_iid_values <- stats::rnorm(length(missing_pixels), 0, iideffect_sd)
      terra::values(iid_ras)[missing_pixels] <- missing_iid_values
    }
    linear_pred <- linear_pred + iid_ras
  } else {
    iid_ras <- NULL
  }

    if(link_function == 'logit') {
    prediction_ras <- 1 / (1 + exp(-1 * linear_pred))
  } else if(link_function == 'log') {
    prediction_ras <- exp(linear_pred)
  } else if(link_function == 'identity') {
    prediction_ras <- linear_pred
  }

    predictions <- list(prediction = prediction_ras,
                      field = field_ras,
                      iid = iid_ras,
                      covariates = cov_contribution)

  return(predictions)
}
