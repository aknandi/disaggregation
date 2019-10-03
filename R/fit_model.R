#' Fit the disaggregation model
#' 
#' \emph{fit_model} function takes a \emph{disag.data} object created by \emph{disaggregation::prepare_data} and performs a Bayesian disaggregation fit
#' 
#' Takes a \emph{disag.data} object created by \emph{disaggregation::prepare_data}
#' 
#' Specify priors for the regression parameters, field and iid effect as a single list.
#' 
#' The \emph{family} and \emph{link} arguments are used to specify the likelihood and link function respectively. 
#' The likelihood function can be one of \emph{gaussian}, \emph{poisson} or \emph{binomial}. 
#' The link function can be one of \emph{logit}, \emph{log} or \emph{identity}.
#' These are specified as strings
#' 
#' The field and iid effect can be turned on or off via the \emph{field} and \emph{iid} logical flags. Default TRUE.
#' 
#' The \emph{its} argument specifies the maximum number of iterations the model can run for to find an optimal point.
#' 
#' The \emph{silent} argument can be used to publish/supress verbose output. Default TRUE.
#' 
#'
#' @param data disag.data object returned by \emph{prepare_data} function that contains all the necessary objects for the model fitting
#' @param priors list of prior values
#' @param family likelihood function: gaussian, binomial or poisson
#' @param link link function: logit, log or identity
#' @param its number of iterations to run the optimisation for
#' @param field logical. Flag the spatial field on or off
#' @param iid logical. Flag the iid effect on or off
#' @param silent logical. Suppress verbose output.
#' 
#' 
#' @return A list is returned of class \code{fit.result}. 
#' The functions \emph{summary}, \emph{print} and \emph{plot} can be used on \code{fit.result}. 
#' The list  of class \code{fit.result} contains:
#'  \item{obj }{The TMB model object returned by \emph{TMB::MakeADFun}.} 
#'  \item{opt }{The optimized model object.} 
#'  \item{sd_out }{The TMB object return by \emph{TMB::sdreport}.}
#'  \item{data }{The \emph{disag.data} object used as an input to the model.}
#'  \item{model_setup }{A list of information on the model setup. Likelihood function (\emph{family}), link function(\emph{link}), logical: whether a field was used (\emph{field}) and logical: whether an iid effect was used (\emph{iid}).}
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

fit_model <- function(data, 
                      priors = NULL, 
                      family = 'gaussian', 
                      link = 'identity', 
                      its = 10, 
                      field = TRUE, 
                      iid = TRUE,
                      silent = TRUE) {
  
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
                         priorsd_log_tau = 2.0)
  
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
    silent = silent,
    DLL = "disaggregation")
  
  message('Fitting model. This may be slow.')
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = its, trace = 0))
  
  
  
  sd_out <- TMB::sdreport(obj, getJointPrecision = TRUE)
  
  if(opt$convergence != 0) warning('The model did not converge. Try increasing its')
  
  model_output <- list(obj = obj,
                       opt = opt,
                       sd_out = sd_out,
                       data = data,
                       model_setup = list(family = family_id, link = link_id, field = field, iid = iid))
  
  class(model_output) <- c('fit.result', 'list')
  
  return(model_output)
  
  
}
