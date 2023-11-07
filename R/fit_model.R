#' Fit the disaggregation model
#'
#' \emph{fit_model} function takes a \emph{disag_data} object created by
#' \code{\link{prepare_data}} and performs a Bayesian disaggregation fit.
#'
#' \strong{The model definition}
#'
#' The disaggregation model makes predictions at the pixel level:
#' \deqn{link(pred_i) = \beta_0 + \beta X + GP(s_i) + u_i}{ link(predi) = \beta 0 + \beta X + GP + u}
#'
#' And then aggregates these predictions to the polygon level using the weighted sum (via the aggregation raster, \eqn{agg_i}{aggi}):
#' \deqn{cases_j = \sum_{i \epsilon j} pred_i \times agg_i}{ casesj = \sum (predi x aggi)}
#' \deqn{rate_j = \frac{\sum_{i \epsilon j} pred_i \times agg_i}{\sum_{i \epsilon j} agg_i}}{ratej = \sum(predi x aggi) / \sum (aggi)}
#'
#' The different likelihood correspond to slightly different models (\eqn{y_j}{yi} is the response count data):
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
#' where \eqn{P(\rho < \rho_{min}) = \rho_{prob}}, and \eqn{\sigma_{min}} and \eqn{\sigma_{prob}} for the variation of the field
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
#' The \emph{silent} argument can be used to publish/suppress verbose output. Default TRUE.
#'
#'
#' @param data disag_data object returned by \code{\link{prepare_data}} function that contains all the necessary objects for the model fitting
#' @param priors list of prior values
#' @param family likelihood function: \emph{gaussian}, \emph{binomial} or \emph{poisson}
#' @param link link function: \emph{logit}, \emph{log} or \emph{identity}
#' @param iterations number of iterations to run the optimisation for
#' @param field logical. Flag the spatial field on or off
#' @param iid logical. Flag the iid effect on or off
#' @param hess_control_parscale Argument to scale parameters during the calculation of the Hessian.
#' Must be the same length as the number of parameters. See \code{\link[stats]{optimHess}} for details.
#' @param hess_control_ndeps Argument to control step sizes during the calculation of the Hessian.
#' Either length 1 (same step size applied to all parameters) or the same length as the number of parameters.
#' Default is 1e-3, try setting a smaller value if you get NaNs in the standard error of the parameters.
#' See \code{\link[stats]{optimHess}} for details.
#' @param silent logical. Suppress verbose output.
#'
#' @return A list is returned of class \code{disag_model}.
#' The functions \emph{summary}, \emph{print} and \emph{plot} can be used on \code{disag_model}.
#' The list  of class \code{disag_model} contains:
#'  \item{obj }{The TMB model object returned by \code{\link[TMB]{MakeADFun}}.}
#'  \item{opt }{The optimized model object returned by \code{\link[stats]{nlminb}}.}
#'  \item{sd_out }{The TMB object returned by \code{\link[TMB]{sdreport}}.}
#'  \item{data }{The \emph{disag_data} object used as an input to the model.}
#'  \item{model_setup }{A list of information on the model setup. Likelihood function (\emph{family}), link function(\emph{link}), logical: whether a field was used (\emph{field}) and logical: whether an iid effect was used (\emph{iid}).}
#'
#' @name fit_model
#' @references Nanda et al. (2023) disaggregation: An R Package for Bayesian
#' Spatial Disaggregation Modeling. <doi:10.18637/jss.v106.i11>
#'
#' @examples
#' \dontrun{
#' polygons <- list()
#' n_polygon_per_side <- 10
#' n_polygons <- n_polygon_per_side * n_polygon_per_side
#' n_pixels_per_side <- n_polygon_per_side * 2
#'
#' for(i in 1:n_polygons) {
#'   row <- ceiling(i/n_polygon_per_side)
#'   col <- ifelse(i %% n_polygon_per_side != 0, i %% n_polygon_per_side, n_polygon_per_side)
#'   xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
#'   polygons[[i]] <- list(cbind(c(xmin, xmax, xmax, xmin, xmin),
#'                               c(ymax, ymax, ymin, ymin, ymax)))
#' }
#'
#' polys <- lapply(polygons,sf::st_polygon)
#' N <- floor(runif(n_polygons, min = 1, max = 100))
#' response_df <- data.frame(area_id = 1:n_polygons, response = runif(n_polygons, min = 0, max = 1000))
#'
#' spdf <- sf::st_sf(response_df, geometry = polys)
#'
#' # Create raster stack
#' r <- terra::rast(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
#' terra::ext(r) <- terra::ext(spdf)
#' r[] <- sapply(1:terra::ncell(r), function(x){
#' rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side), 3))}
#' r2 <- terra::rast(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
#' terra::ext(r2) <- terra::ext(spdf)
#' r2[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ceiling(x/n_pixels_per_side), 3))
#' cov_stack <- c(r, r2)
#' names(cov_stack) <- c('layer1', 'layer2')
#'
#' test_data <- prepare_data(polygon_shapefile = spdf,
#'                           covariate_rasters = cov_stack)
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
                      hess_control_parscale = NULL,
                      hess_control_ndeps = 1e-4,
                      silent = TRUE) {

  .Deprecated(new = 'disag_model', msg = "'fit_model' will be removed in the next version. Please use 'disag_model' instead")

  model_output <- disag_model(data,
                              priors = priors,
                              family = family,
                              link = link,
                              iterations = iterations,
                              field = field,
                              iid = iid,
                              hess_control_parscale = hess_control_parscale,
                              hess_control_ndeps = hess_control_ndeps,
                              silent = silent)

  return(model_output)


}

#' @export
#' @rdname fit_model

disag_model <- function(data,
                        priors = NULL,
                        family = 'gaussian',
                        link = 'identity',
                        iterations = 100,
                        field = TRUE,
                        iid = TRUE,
                        hess_control_parscale = NULL,
                        hess_control_ndeps = 1e-4,
                        silent = TRUE) {


  stopifnot(inherits(data, 'disag_data'))
  if(!is.null(priors)) stopifnot(inherits(priors, 'list'))
  stopifnot(inherits(iterations, 'numeric'))

  obj <- make_model_object(data = data,
                           priors = priors,
                           family = family,
                           link = link,
                           field = field,
                           iid = iid,
                           silent = silent)

  message('Fitting model. This may be slow.')
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = iterations, trace = 0))

  if(opt$convergence != 0) warning('The model did not converge. Try increasing the number of iterations')

  # Get hess control parameters into a list.
  hess_control <- setup_hess_control(opt, hess_control_parscale, hess_control_ndeps)

  # Calculate the hessian
  hess <- stats::optimHess(opt$par, fn = obj$fn, gr = obj$gr, control = hess_control)

  # Calc uncertainty using the fixed hessian from above.
  sd_out <- TMB::sdreport(obj, getJointPrecision = TRUE, hessian.fixed = hess)


  sd_out <- TMB::sdreport(obj, getJointPrecision = TRUE)

  # Rename parameters to match layers
  # Need to change in sd_out as well
  # names(opt$par)[names(opt$par) == 'slope'] <- names(data$covariate_rasters)

  model_output <- list(obj = obj,
                       opt = opt,
                       sd_out = sd_out,
                       data = data,
                       model_setup = list(family = family, link = link, field = field, iid = iid))

  class(model_output) <- c('disag_model', 'list')

  return(model_output)
}

#' Create the TMB model object for the disaggregation model
#'
#' \emph{make_model_object} function takes a \emph{disag_data} object created by \code{\link{prepare_data}}
#' and creates a TMB model object to be used in fitting.
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
#' The different likelihood correspond to slightly different models (\eqn{y_j}{yi} is the response count data):
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
#' Specify priors for the regression parameters, field and iid effect as a single named list. Hyperpriors for the field
#' are given as penalised complexity priors you specify \eqn{\rho_{min}} and \eqn{\rho_{prob}} for the range of the field
#' where \eqn{P(\rho < \rho_{min}) = \rho_{prob}}, and \eqn{\sigma_{min}} and \eqn{\sigma_{prob}} for the variation of the field
#' where \eqn{P(\sigma > \sigma_{min}) = \sigma_{prob}}. Also, specify pc priors for the iid effect.
#'
#' The precise names and default values for these priors are:
#' \itemize{
#' \item priormean_intercept: 0
#' \item priorsd_intercept: 10.0
#' \item priormean_slope: 0.0
#' \item priorsd_slope: 0.5
#' \item prior_rho_min: A third the length of the diagonal of the bounding box.
#' \item prior_rho_prob: 0.1
#' \item prior_sigma_max: sd(response/mean(response))
#' \item prior_sigma_prob: 0.1
#' \item prior_iideffect_sd_max: 0.1
#' \item prior_iideffect_sd_prob: 0.01
#' }
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
#' @param data disag_data object returned by \code{\link{prepare_data}} function that contains all the necessary objects for the model fitting
#' @param priors list of prior values
#' @param family likelihood function: \emph{gaussian}, \emph{binomial} or \emph{poisson}
#' @param link link function: \emph{logit}, \emph{log} or \emph{identity}
#' @param field logical. Flag the spatial field on or off
#' @param iid logical. Flag the iid effect on or off
#' @param silent logical. Suppress verbose output.
#'
#' @return The TMB model object returned by \code{\link[TMB]{MakeADFun}}.
#'
#' @name make_model_object
#'
#' @examples
#' \dontrun{
#' polygons <- list()
#' n_polygon_per_side <- 10
#' n_polygons <- n_polygon_per_side * n_polygon_per_side
#' n_pixels_per_side <- n_polygon_per_side * 2
#'
#' for(i in 1:n_polygons) {
#'   row <- ceiling(i/n_polygon_per_side)
#'   col <- ifelse(i %% n_polygon_per_side != 0, i %% n_polygon_per_side, n_polygon_per_side)
#'   xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
#'   polygons[[i]] <- list(cbind(c(xmin, xmax, xmax, xmin, xmin),
#'                               c(ymax, ymax, ymin, ymin, ymax)))
#' }
#'
#' polys <- lapply(polygons,sf::st_polygon)
#' N <- floor(runif(n_polygons, min = 1, max = 100))
#' response_df <- data.frame(area_id = 1:n_polygons, response = runif(n_polygons, min = 0, max = 1000))
#'
#' spdf <- sf::st_sf(response_df, geometry = polys)
#'
#' # Create raster stack
#' r <- terra::rast(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
#' terra::ext(r) <- terra::ext(spdf)
#' r[] <- sapply(1:terra::ncell(r), function(x){
#' rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side), 3))}
#' r2 <- terra::rast(ncol=n_pixels_per_side, nrow=n_pixels_per_side)
#' terra::ext(r2) <- terra::ext(spdf)
#' r2[] <- sapply(1:terra::ncell(r), function(x) rnorm(1, ceiling(x/n_pixels_per_side), 3))
#' cov_stack <- c(r, r2)
#' names(cov_stack) <- c('layer1', 'layer2')
#'
#' test_data <- prepare_data(polygon_shapefile = spdf,
#'                           covariate_rasters = cov_stack)
#'
#'  result <- make_model_object(test_data)
#'  }
#'
#' @export
#'

make_model_object <- function(data,
                              priors = NULL,
                              family = 'gaussian',
                              link = 'identity',
                              field = TRUE,
                              iid = TRUE,
                              silent = TRUE) {


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
    warning('You are using both a gaussian likelihood and an iid effect. Using both of these is redundant as they are
            having the same effect on the model. Consider setting iid = FALSE.')
  }

  if(is.null(data$mesh)) {
    stop('Your data object must contain an INLA mesh.')
  }

  nu = 1
  # Sort out mesh bits
  spde <- (INLA::inla.spde2.matern(data$mesh, alpha = nu + 1)$param.inla)[c("M0", "M1", "M2")]
  Apix <- fmesher::fm_evaluator(data$mesh, loc = data$coordsForFit)$proj$A
  n_s <- nrow(spde$M0)

  cov_matrix <- as.matrix(data$covariate_data[, (names(data$covariate_data) %in% names(data$covariate_rasters))])
  # If we have exactly one column we don't have to transpose. Sure this
  #   this could be cleaner but I don't know how.
  if(ncol(cov_matrix) == 1){
    cov_matrix <- as.matrix(apply(cov_matrix, 1, as.numeric))
  } else {
    cov_matrix <- t(apply(cov_matrix, 1, as.numeric))
  }

  # Construct sensible default field hyperpriors
  limits <- sf::st_bbox(data$polygon_shapefile)
  hypotenuse <- sqrt((limits$xmax - limits$xmin)^2 + (limits$ymax - limits$ymin)^2)
  prior_rho <- hypotenuse/3

  prior_sigma <- sd(data$polygon_data$response/mean(data$polygon_data$response))

  # Default priors if they are not specified
  default_priors <- list(priormean_intercept = 0,
                         priorsd_intercept = 10.0,
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
                     nu = nu,
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

  return(obj)
}



# Setup hessian control
setup_hess_control <- function(opt,hess_control_parscale, hess_control_ndeps){
  hess_control <- list()
  # hess_control_parscale should always either be null or a vector of the correct length.
  if(!is.null(hess_control_parscale)){
    if(length(hess_control_parscale) != length(opt$par)){
      stop(paste('hess_control_parscale must either be NULL or a vector of length', length(opt$par)))
    }
    hess_control$parscale <- hess_control_parscale
  }
  # hess_control_ndeps can either be length 1 (default) or correct length vecot.
  if(length(hess_control_ndeps) == 1){
    hess_control$ndeps <- rep(hess_control_ndeps, length(opt$par))
  } else {
    if(length(hess_control_ndeps) != length(opt$par)){
      stop(paste('hess_control_ndeps must either be NULL or a vector of length', length(opt$par)))
    }
    hess_control$ndeps <- hess_control_ndeps
  }
  return(hess_control)
}
