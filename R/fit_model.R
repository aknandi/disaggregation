#' Fit the disaggregation model
#' 
#' \emph{fit_model} function takes a \emph{disag.data} object created by \code{\link{prepare_data}} and performs a Bayesian disaggregation fit.
#' 
#' \strong{The model definition}
#' 
#' The disaggregation model make predictions at the pixel level:
#' \deqn{link(pred_i) = \beta_0 + \beta X + GP(s_i) + u_i}{ link(predi) = \beta 0 + \beta X + GP + u}
#' 
#' And then aggregates these predictions to the polygon level using the weighted sum (via the aggregation raster, \eqn{agg_i}{aggi}):
#' \deqn{cases_j = \sum_{i \epsilon j} pred_i \times agg_i}{ casesj = \sum (predi x aggi)}
#' \deqn{rate_j = \frac{\sum_{i \epsilon j} pred_i \times agg_i}{\sum_{i \epsilon j} agg_i}}{ratej = \sum(predi x aggi) / \sum (aggi)}
#' 
#' The different likelihood correspond to slightly different models (\eqn{y_j}{yi} is the repsonse count data):
#' \itemize{
#'   \item Gaussian: 
#'    If \eqn{\sigma} is the dispersion of the pixel data, \eqn{\sigma_j}{\sigmaj} is the dispersion of the polygon data, where 
#'    \eqn{\sigma_j = \sigma \sqrt{\sum agg_i^2} / \sum agg_i }{\sigmaj = \sigma x { \sqrt \sum (aggi ^ 2) } / \sum aggi}
#'    \deqn{dnorm(y_j/\sum agg_i, rate_j, \sigma_j)}{dnorm(yj / \sum aggi, ratej, \sigmaj)} - predicts incidence rate.
#'   \item Binomial: 
#'    For a survey in polygon j, \eqn{y_j}{yj} is the number positive and \eqn{N_j}{Nj} is the number tested.
#'    \deqn{dbinom(y_j, N_j, rate_j)}{dbinom(yj, Nj, ratej)} - predicts prevalence rate.
#'   \item Poisson: 
#'    \deqn{dpois(y_j, cases_j)}{dpois(yj, casesj)} - predicts incidence count.
#' }
#' 
#' Specify priors for the regression parameters, field and iid effect as a single list. Hyperpriors for the field 
#' are given as penalised complexity priors you specify \eqn{\rho_{min}} and \eqn{\rho_{prob}} for the range of the field 
#' where \eqn{P(\rho < \rho_{min}) = \rho_{prob}}, and \eqn{\sigma_{min}$ and $\sigma_{prob}} for the variation of the field 
#' where \eqn{P(\sigma > \sigma_{min}) = \sigma_{prob}}. Also, specify pc priors for the iid effect
#' 
#' The \emph{family} and \emph{link} arguments are used to specify the likelihood and link function respectively. 
#' The likelihood function can be one of \emph{gaussian}, \emph{poisson} or \emph{binomial}. 
#' The link function can be one of \emph{logit}, \emph{log} or \emph{identity}.
#' These are specified as strings.
#' 
#' The field and iid effect can be turned on or off via the \emph{field} and \emph{iid} logical flags. Both are default TRUE.
#' 
#' The \emph{iterations} argument specifies the maximum number of iterations the model can run for to find an optimal point.
#' 
#' The \emph{silent} argument can be used to publish/supress verbose output. Default TRUE.
#' 
#'
#' @param data disag.data object returned by \code{\link{prepare_data}} function that contains all the necessary objects for the model fitting
#' @param priors list of prior values
#' @param family likelihood function: \emph{gaussian}, \emph{binomial} or \emph{poisson}
#' @param link link function: \emph{logit}, \emph{log} or \emph{identity}
#' @param iterations number of iterations to run the optimisation for
#' @param field logical. Flag the spatial field on or off
#' @param iid logical. Flag the iid effect on or off
#' @param silent logical. Suppress verbose output.
#' 
#' 
#' @return A list is returned of class \code{fit.result}. 
#' The functions \emph{summary}, \emph{print} and \emph{plot} can be used on \code{fit.result}. 
#' The list  of class \code{fit.result} contains:
#'  \item{obj }{The TMB model object returned by \code{\link[TMB]{MakeADFun}}.} 
#'  \item{opt }{The optimized model object returned by \code{\link[stats]{nlminb}}.} 
#'  \item{sd_out }{The TMB object returned by \code{\link[TMB]{sdreport}}.}
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
#'  result <- fit_model(test_data, iterations = 2)
#'  }
#' 
#' @export

fit_model <- function(data, 
                      priors = NULL, 
                      family = 'gaussian', 
                      link = 'identity', 
                      iterations = 100, 
                      field = TRUE, 
                      iid = TRUE,
                      silent = TRUE) {
  
  stopifnot(inherits(data, 'disag.data'))
  if(!is.null(priors)) stopifnot(inherits(priors, 'list'))
  stopifnot(inherits(iterations, 'numeric'))
  
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
  
  if(family == 'gaussian' & iid) {
    warning('You are using both a gaussian likeihood and an iid effect. Using both of these is redundant as they are 
            having the same effect on the model. Consider setting iid = FALSE.')
  }
  
  if(is.null(data$mesh)) {
    stop('Your data object must contain an INLA mesh.')
  }
  
  # Sort out mesh bits
  spde <- (INLA::inla.spde2.matern(data$mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  Apix <- INLA::inla.mesh.project(data$mesh, loc = data$coordsForFit)$A
  n_s <- nrow(spde$M0)
  
  cov_matrix <- as.matrix(data$covariate_data[, -c(1:2)])
  cov_matrix <- t(apply(cov_matrix, 1,as.numeric))
  
  # Construct sensible default field hyperpriors
  limits <- sp::bbox(data$polygon_shapefile)
  hypontenuse <- sqrt((limits[1,2] - limits[1,1])^2 + (limits[2,2] - limits[2,1])^2)
  prior_rho <- hypontenuse/3
  
  prior_sigma <- sd(data$polygon_data$response/mean(data$polygon_data$response))
  
  # Default priors if they are not specified
  default_priors <- list(priormean_intercept = -4.0,
                         priorsd_intercept = 2.0,
                         priormean_slope = 0.0,
                         priorsd_slope = 0.5,
                         prior_rho_min = prior_rho,
                         prior_rho_prob = 0.1,
                         prior_sigma_max = prior_sigma,
                         prior_sigma_prob = 0.1,
                         prior_iideffect_sd_max = 0.1,
                         prior_iideffect_sd_prob = 0.01)
  
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
                     log_tau_gaussian = 8,
                     iideffect = rep(0, nrow(data$polygon_data)),
                     iideffect_log_tau = 1,
                     log_sigma = 0,
                     log_rho = 4,
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
    tmb_map <- c(tmb_map, list(log_sigma = as.factor(NA),
                               log_rho = as.factor(NA),
                               nodemean = factor(rep(NA, n_s))))
  }
  if(!iid) {
    tmb_map <- c(tmb_map, list(iideffect_log_tau = as.factor(NA),
                               iideffect = factor(rep(NA, nrow(data$polygon_data)))))
  }
  if(family_id != 0) { # if not gaussian do not need a dispersion in likelihood
    tmb_map <- c(tmb_map, list(log_tau_gaussian = as.factor(NA)))
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
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = iterations, trace = 0))
  
  
  
  sd_out <- TMB::sdreport(obj, getJointPrecision = TRUE)
  
  if(opt$convergence != 0) warning('The model did not converge. Try increasing the number of iterations')
  
  model_output <- list(obj = obj,
                       opt = opt,
                       sd_out = sd_out,
                       data = data,
                       model_setup = list(family = family, link = link, field = field, iid = iid))
  
  class(model_output) <- c('fit.result', 'list')
  
  return(model_output)
  
  
}
