#' Function to fit the disaggregation model
#'
#' @param data disag.data object returned by prepare_data function that contains all the necessary objects for the model fitting
#' @param priors list of prior values
#' @param family likelihood function: gaussian, binomial or poisson
#' @param link link function: logit, log or identity
#' @param its number of iterations to run the optimisation for
#' @param field boolean to flag spatial field on or off
#' @param iid boolean to flag iid effect on or off
#' 
#' @name fit_model
#'
#' @examples 
#' \dontrun{
#'  polygons <- list()
#'  for(i in 1:100) {
#'   row <- ceiling(i/10)
#'   col <- ifelse(i %% 10 != 0, i %% 10, 10)
#'   xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
#'   polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
#'  }
#' 
#'  polys <- do.call(raster::spPolygons, polygons)
#'  response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
#'  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
#' 
#'  r <- raster::raster(ncol=20, nrow=20)
#'  r <- raster::setExtent(r, raster::extent(spdf))
#'  r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
#'  r2 <- raster::raster(ncol=20, nrow=20)
#'  r2 <- raster::setExtent(r2, raster::extent(spdf))
#'  r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
#'  cov_rasters <- raster::stack(r, r2)
#' 
#'  cl <- parallel::makeCluster(2)
#'  doParallel::registerDoParallel(cl)
#'  test_data <- prepare_data(polygon_shapefile = spdf, 
#'                            covariate_rasters = cov_rasters)
#'  parallel::stopCluster(cl)
#'  foreach::registerDoSEQ()
#'                          
#'  result <- fit_model(test_data, its = 2)
#'  }
#' 
#' @export

fit_model <- function(data, priors = NULL, family = 'gaussian', link = 'identity', its = 10, field = TRUE, iid = TRUE) {
  
  stopifnot(inherits(data, 'disag.data'))
  if(!is.null(priors)) stopifnot(inherits(priors, 'list'))
  stopifnot(inherits(its, 'numeric'))
  
  # Check that binomial model has sample_size values supplied
  if(family == 'binomial') {
    if(sum(is.na(data$polygon_data$N)) != 0) {
      stop("There are NAs in the sample sizes. These must be supplied for a binomial likelihood")
    }
  }
  
  if(family == 'gaussian') {
    family_id = 0
  } else if(family == 'binomial') {
    family_id = 1
  } else if(family == 'poisson') {
    family_id = 2
  } else {
    stop(paste(family, "is not a valid likelihood"))
  }
  
  if(link == 'logit') {
    link_id = 0
  } else if(link == 'log') {
    link_id = 1
  } else if(link == 'identity') {
    link_id = 2
  } else {
    stop(paste(link, "is not a valid link function"))
  }
  
  # Sort out mesh bits
  spde <- (INLA::inla.spde2.matern(data$mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  Apix <- INLA::inla.mesh.project(data$mesh, loc = data$coords)$A
  n_s <- nrow(spde$M0)
  
  cov_matrix <- as.matrix(data$covariate_data[, -c(1:2)])
  cov_matrix <- t(apply(cov_matrix, 1,as.numeric))
  
  # Default priors if they are not specified
  default_priors <- list(polygon_sd_mean = 0.1,
                         polygon_sd_sd = 0.1,
                         priormean_intercept = -4.0,
                         priorsd_intercept = 2.0,
                         priormean_slope = 0.0,
                         priorsd_slope = 0.5,
                         priorsd_iideffect = 0.05,
                         priormean_log_kappa = -3,
                         priorsd_log_kappa = 0.5,
                         priormean_log_tau = -0.50,
                         priorsd_log_tau = 2.0,
                         priorsd_log_tau = 3.0)
  
  # Replace with any specified priors
  if(!is.null(priors)) {
    final_priors <- default_priors
    prior_names <- names(priors)
    # Check all input priors are named correctly
    if(sum(!(prior_names %in% names(final_priors))) != 0) {
      stop(paste(prior_names[!(prior_names %in% names(final_priors))], 'is not the name of a prior'))
    }
    # Check priors are not specified multiple times
    if(sum(duplicated(names(priors))) != 0) {
      message(paste(names(priors)[duplicated(names(priors))],'are specified multiple times. Will only take first value'))
      priors <- priors[!duplicated(names(priors))]
    }
    # Replace default value with new prior value
    for(i in 1:length(priors)) {
      prior_to_replace <- prior_names[i]
      final_priors[[prior_to_replace]] <- priors[[prior_to_replace]]
    }
  } else {
    final_priors <- default_priors
  }

  parameters <- list(intercept = -5,
                     slope = rep(0, ncol(cov_matrix)),
                     polygon_sd = 0.1,
                     iideffect = rep(0, nrow(data$polygon_data)),
                     log_kappa = -3,
                     log_tau = -0.5,
                     nodemean = rep(0, n_s))
  
  input_data <- list(x = cov_matrix,
                     aggregation_values = data$aggregation_pixels,
                     Apixel = Apix,
                     spde = spde,
                     startendindex = data$startendindex,
                     polygon_response_data = data$polygon_data$response,
                     response_sample_size = data$polygon_data$N,
                     family = family_id,
                     link = link_id,
                     field = as.integer(field),
                     iid = as.integer(iid))
  
  input_data <- c(input_data, final_priors)
  
  tmb_map <- list()
  if(!field) {
    tmb_map <- c(tmb_map, list(log_kappa = as.factor(NA),
                               log_tau = as.factor(NA),
                               nodemean = factor(rep(NA, n_s))))
  }
  if(!iid) {
    tmb_map <- c(tmb_map, list(iideffect = factor(rep(NA, nrow(data$polygon_data)))))
  }
  if(family_id != 0) { # if not gaussian do not need a dispersion in likelihood
    tmb_map <- c(tmb_map, list(polygon_sd = as.factor(NA)))
  }
  
  random_effects <- c()
  if(field) {
    random_effects <- c(random_effects, 'nodemean')
  }
  if(iid) {
    random_effects <- c(random_effects, 'iideffect')
  }
  
  obj <- TMB::MakeADFun(
    data = input_data, 
    parameters = parameters,
    map = tmb_map,
    random = random_effects,
    DLL = "disaggregation")
  
  
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = its, eval.max = 2*its, trace = 0))
  
  sd_out <- TMB::sdreport(obj, getJointPrecision = TRUE)
  
  model_output <- list(obj = obj,
                       opt = opt,
                       sd_out = sd_out,
                       data = data,
                       model_setup = list(family = family_id, link = link_id, field = field, iid = iid))
  
  class(model_output) <- c('fit.result', 'list')
  
  return(model_output)
  
  
}
